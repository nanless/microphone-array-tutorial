"""Measure local UDP delivery of a LibriWASN excerpt and read a public DUST trace."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import io
from itertools import islice
import json
from pathlib import Path
import socket
import struct
import time

import numpy as np

from codes.chapters.ch13.examples.clock_recovery import (
    BLOCK, FS, GAP, _libriwasn_inputs, _run_pair,
)


PACKET = 256
PACKET_NS = PACKET * 1_000_000_000 // FS
GRACE_NS = 2 * PACKET_NS
GAP_PACKETS = range(GAP[0] * BLOCK // PACKET, GAP[1] * BLOCK // PACKET)
LATE_PACKET = GAP_PACKETS.stop - 1
HEADER = struct.Struct("!IIQ")
DUST_ROWS = 20_000
DUST_BYTES = 24_265_783
DUST_MD5 = "957e72c6e89b570ec68586e71b7e1272"
DUST_PREFIX_SHA256 = "56c13a7a4054e1de8be9d20cb752f8ca345fff703517ad3278a55d5e22a82fad"


def _packet_time(zero, packet):
    return zero + (packet + 1) * PACKET_NS


def _send(sender, peer, zero, pcm):
    """Use media sample indices even when suppressing or delaying a datagram."""
    count = len(pcm) // PACKET
    late_at = _packet_time(zero, LATE_PACKET) + 4 * PACKET_NS
    schedule = [(_packet_time(zero, n), n) for n in range(count) if n not in GAP_PACKETS]
    schedule.append((late_at, LATE_PACKET))
    for due, packet in sorted(schedule):
        while time.monotonic_ns() < due:
            time.sleep(max(0, due - time.monotonic_ns()) / 1e9)
        start = packet * PACKET
        body = pcm[start:start + PACKET].tobytes()
        sender.sendto(HEADER.pack(packet, start, time.monotonic_ns()) + body, peer)


def _receive(receiver, zero, count):
    pcm = np.zeros(count * PACKET, dtype="<i2")
    valid = np.zeros(pcm.size, dtype=bool)
    arrivals = {}
    end = _packet_time(zero, count - 1) + 250_000_000
    while time.monotonic_ns() < end:
        receiver.settimeout(max(0.001, (end - time.monotonic_ns()) / 1e9))
        try:
            data = receiver.recv(HEADER.size + 2 * PACKET + 1)
        except TimeoutError:
            break
        received = time.monotonic_ns()
        if len(data) != HEADER.size + 2 * PACKET:
            raise ValueError("unexpected UDP datagram length")
        number, start, sent = HEADER.unpack_from(data)
        if number >= count or start != number * PACKET or number in arrivals or sent > received:
            raise ValueError("unexpected packet number, sample index, or send time")
        deadline = _packet_time(zero, number) + GRACE_NS
        on_time = received <= deadline
        arrivals[number] = {"first_sample": start, "send_ns": sent - zero, "arrival_ns": received - zero,
                            "deadline_ns": deadline - zero, "accepted": on_time}
        if on_time:
            valid[start:start + PACKET] = True
            pcm[start:start + PACKET] = np.frombuffer(data, dtype="<i2", count=PACKET, offset=HEADER.size)
    return pcm, valid, arrivals


def _local_udp(recording):
    samples, metadata = _libriwasn_inputs(recording)
    device = "Pixel6a"
    original = (samples[device] * 32768).astype("<i2")
    expected = original.size // PACKET
    with (socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as sender,
          socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as receiver,
          ThreadPoolExecutor(max_workers=1) as worker):
        receiver.bind(("127.0.0.1", 0))
        receiver.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, 256 * 1024)
        zero = time.monotonic_ns() + 100_000_000
        sending = worker.submit(_send, sender, receiver.getsockname(), zero, original)
        pcm, accepted, arrivals = _receive(receiver, zero, expected)
        sending.result()
    if not np.array_equal(original[accepted], pcm[accepted]):
        raise ValueError("received on-time samples differ from the original PCM")
    missing = sorted(set(range(expected)) - arrivals.keys())
    late = sorted(n for n, row in arrivals.items() if not row["accepted"])
    sample_arrival_ns = np.zeros(pcm.size, dtype=np.int64)
    for row in arrivals.values():
        start = row["first_sample"]
        sample_arrival_ns[start:start + PACKET] = row["arrival_ns"]
    timing = _run_pair(samples[metadata["reference"]], pcm.astype(float) / 32768,
                       received=accepted, arrival_ns=sample_arrival_ns)
    return {
        "recording": metadata, "remote": device, "network": "local IPv4 UDP loopback",
        "packet_samples": PACKET, "packet_period_ns": PACKET_NS,
        "receiver_grace_ns_from_nominal_packet_end": GRACE_NS,
        "expected_packets": expected, "accepted_on_time": expected - len(missing) - len(late),
        "omitted_by_sender": list(range(GAP_PACKETS.start, LATE_PACKET)),
        "delayed_by_sender": LATE_PACKET, "missing_at_receiver": missing, "received_late": late,
        "arrivals": arrivals, "alignment_replay": timing,
    }


def _dust(path):
    full_file = path.stat().st_size == DUST_BYTES
    if full_file:
        with path.open("rb") as stream:
            if hashlib.file_digest(stream, "md5").hexdigest() != DUST_MD5:
                raise ValueError("DUST archive does not match the published MD5")
    # The analyzed prefix is independently fixed, including when only the beginning of
    # the gzip download is available. Text mode normalizes line endings to LF.
    with gzip.open(path, "rt", encoding="ascii") as stream:
        excerpt = "".join(islice(stream, DUST_ROWS + 1))
    digest = hashlib.sha256(excerpt.encode("ascii")).hexdigest()
    if digest != DUST_PREFIX_SHA256:
        raise ValueError("expected the fixed prefix of Zenodo sanri_14A.csv.gz")
    reader = csv.DictReader(io.StringIO(excerpt))
    fields = ["read_index", "send_index", "tstamp_ns"]
    rows = [tuple(int(row[field]) for field in fields) for row in reader]
    step = 24
    if len(rows) != DUST_ROWS or any(sent != i * step for i, (_, sent, _) in enumerate(rows)):
        raise ValueError("the selected DUST rows no longer have consecutive 24-sample packets")
    # DUST paper, Eq. (2): W = B + S / (1 + d / 1e6). Here d = 0 in both controls.
    def replay(buffer_samples):
        write = [buffer_samples + sent for _, sent, _ in rows]
        missed = [min(step, max(0, read - start)) for (read, _, _), start in zip(rows, write)]
        intact, fully = missed.count(0), missed.count(step)
        return missed, {
            "write_offset_samples": buffer_samples, "drift_correction_ppm": 0,
            "intact_arrivals": intact, "partially_late_arrivals": len(rows) - intact - fully,
            "fully_late_arrivals": fully, "total_samples_unusable": sum(missed),
        }

    missed, no_buffer = replay(0)
    _, one_packet_buffer = replay(step)
    first_lost = missed.index(step)
    first_usable = next(i for i in range(first_lost + 1, len(missed)) if missed[i] < step)
    first_intact = next(i for i in range(first_usable, len(missed)) if missed[i] == 0)

    def packet(i):
        return dict(zip(fields, rows[i])) | {"zero_based_row": i, "samples_unusable": missed[i]}

    return {
        "doi": "10.5281/zenodo.12726671", "member": "sanri_14A.csv.gz", "published_md5": DUST_MD5,
        "full_archive_md5_verified": full_file, "selected_csv_prefix_sha256": digest,
        "sample_rate_hz": 44_100, "packet_samples_in_selected_rows": step, "rows": len(rows),
        "skipped_sender_packet_indices": 0, "playback_replays": [no_buffer, one_packet_buffer],
        "no_buffer_example": [packet(first_lost), packet(first_usable), packet(first_intact)],
        "first_intact_elapsed_ns": rows[first_intact][2] - rows[first_lost][2],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--libriwasn", type=Path, help="official ZIP or extracted session8 record directory")
    parser.add_argument("--dust", type=Path, help="the official DUST sanri_14A.csv.gz")
    args = parser.parse_args()
    if args.libriwasn is None and args.dust is None:
        parser.error("give --libriwasn, --dust, or both")
    result = {}
    if args.dust is not None:
        result["dust"] = _dust(args.dust)
    if args.libriwasn is not None:
        result["local_udp"] = _local_udp(args.libriwasn)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
