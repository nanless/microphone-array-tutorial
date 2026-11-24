"""Run or score short, shared AMI far-field meeting excerpts offline."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import dataclass
from hashlib import file_digest, sha256
import gzip
from importlib import import_module
import json
import math
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from time import perf_counter
from types import ModuleType
import wave


SAMPLE_RATE = 16_000
CHANNELS = (1, 2)
WINDOWS = {3: (1_040_000, 1_326_240), 4: (1_327_680, 1_489_600), 21: (4_495_680, 4_818_880)}
CUTSET_SHA256 = "6c6f4d9007df1de49c8ce9fd3072a67b61595608fc043caa5481d32e195f95a4"
PCM_SHA256 = {
    ("ES2011a-3", 1): "ec6a0d335204666236626d6677e60ad3c6ba4252e321b6a085ef3eb0b487ca29",
    ("ES2011a-3", 2): "09729706a5f945b12f9e86e8da556e7b86f17cab3fd5fb32389979435576455f",
    ("ES2011a-4", 1): "7a9d752bb26d77ee31c863cc34c1fd8c1dcc2db8634e1a6a8df4bb41da6005d3",
    ("ES2011a-4", 2): "af2cbfe78c9309b09858189d8b0d92f56bc629567aba5d38619ef6af6521877b",
    ("ES2011a-21", 1): "1321741e9af012add4a6e4293e985c0314a60dfe37058b9192622d1d6eb1ef90",
    ("ES2011a-21", 2): "325686ee72d70291bc6d6eb5c53f58826758edaf502415226f1203a98b146f3d",
}
MLX_WHISPER_SHA256 = "0e03a5993d6eea43b07ee2dcc772b0e4cef5bb227257dacc24bf289387d49186"
MLX_WHISPER_CONFIG_SHA256 = "aaff20ce8f69beddee3fe0cc1e08f4e92f58586cb9f12ba00a6f73cbfec1cb1c"
DICOW_SHA256 = "73dee45746bb4a6e5e99d8553f582a70b21b6f1499dea68433d6a7178eb4183b"
DICOW_SOURCE_SHA256 = {
    "FDDT.py": "fe8b1dec493fdae366b5474de3d1ac0004adcac80036533dacd1ea5f35e0c945",
    "config.py": "fc9eea09841fb6558e68dda42c327332005b5145b144d2fddfe463bff3e9c971",
    "decoding.py": "b713399e01804690b84eea5e4cd58d5c475d8da4217835478a364bff15243f46",
    "encoder.py": "ec63648d28e0950ed4e4b0651bf651c60aab613a3bc603b6a9a97f704dc325ef",
    "generation.py": "a91f671083f385f9dd4469e1f2c24e85fe2ddbd089c572905e5f8259a6cdd38a",
    "layers.py": "e91a8af1e4d62e5b112032e42b89796b3430a0c05c49f8c5031dff6e37c8a842",
    "modeling_dicow.py": "aa335c71078cfad5b84b5e323dd56f232289552a0052b02693819a591e1049b1",
    "utils.py": "f25f0617517f7ed604361643dc875fc5d2d3500376134dbbe240c82688b6714b",
}


@dataclass(frozen=True)
class Clip:
    name: str
    start: int
    end: int
    turns: tuple[tuple[int, int, str, str], ...]

    @property
    def samples(self):
        return self.end - self.start

    @property
    def speakers(self):
        return sorted({speaker for _, _, speaker, _ in self.turns})


def require_sha256(path, expected):
    with path.open("rb") as source:
        found = file_digest(source, "sha256").hexdigest()
    if found != expected:
        raise ValueError(f"SHA-256 mismatch: {path}")


def load_clips(cutset):
    require_sha256(cutset, CUTSET_SHA256)
    expected = {f"sdm_ES2011a-2-0.001-{index}": index for index in WINDOWS}
    clips = {}
    with gzip.open(cutset, "rt", encoding="utf-8") as source:
        for line in source:
            row = json.loads(line)
            if row["id"] not in expected:
                continue
            index = expected[row["id"]]
            start = round(SAMPLE_RATE * row["start"])
            end = start + round(SAMPLE_RATE * row["duration"])
            if (start, end) != WINDOWS[index] or index in clips:
                raise ValueError(f"Unexpected AMI excerpt: {row['id']}")
            turns = []
            for item in row["supervisions"]:
                first = round(SAMPLE_RATE * item["start"])
                last = round(SAMPLE_RATE * (item["start"] + item["duration"]))
                if item["recording_id"] != "ES2011a" or not 0 <= first < last <= end - start:
                    raise ValueError(f"Unexpected AMI supervision: {item['id']}")
                turns.append((first, last, item["speaker"], item["text"]))
            clip = Clip(f"ES2011a-{index}", start, end, tuple(sorted(turns)))
            if len(clip.speakers) != 2:
                raise ValueError(f"Expected two speakers in {clip.name}")
            clips[index] = clip
    if set(clips) != set(WINDOWS):
        raise ValueError("Missing an AMI excerpt")
    return [clips[index] for index in WINDOWS]


def clip_path(directory, clip, channel):
    return directory / f"{clip.name}.Array1-{channel:02}.wav"


def read_pcm(path, start, count, *, exact=False):
    with wave.open(str(path), "rb") as source:
        if (source.getnchannels(), source.getsampwidth(), source.getframerate()) != (1, 2, SAMPLE_RATE):
            raise ValueError(f"Expected 16 kHz mono PCM16: {path}")
        total = source.getnframes()
        if total < start + count or (exact and total != count):
            raise ValueError(f"Incorrect sample count: {path}")
        source.setpos(start)
        data = source.readframes(count)
    if len(data) != 2 * count:
        raise ValueError(f"Truncated PCM data: {path}")
    return data


def require_clip_pcm(data, clip, channel):
    if sha256(data).hexdigest() != PCM_SHA256[clip.name, channel]:
        raise ValueError(f"Unexpected PCM samples: {clip.name}, channel {channel}")
    return data


def prepare(clips, original, destination):
    destination.mkdir(exist_ok=True)
    cuts = {}
    for clip in clips:
        for channel in CHANNELS:
            path = original / f"ES2011a.Array1-{channel:02}.wav"
            data = read_pcm(path, clip.start, clip.samples)
            cuts[clip_path(destination, clip, channel)] = require_clip_pcm(data, clip, channel)
    if any(path.exists() for path in cuts):
        raise FileExistsError("Use a fresh destination for these AMI excerpts")
    for path, data in cuts.items():
        with path.open("xb") as stream, wave.open(stream, "wb") as output:
            output.setparams((1, 2, SAMPLE_RATE, len(data) // 2, "NONE", "not compressed"))
            output.writeframes(data)


def load_audio(clip, directory, channels):
    import numpy as np

    waveforms = []
    for channel in channels:
        data = read_pcm(clip_path(directory, clip, channel), 0, clip.samples, exact=True)
        waveforms.append(np.frombuffer(require_clip_pcm(data, clip, channel), dtype="<i2"))
    return np.stack(waveforms).astype(np.float32) / 32768


def activity(clip, centers):
    import numpy as np

    speakers = clip.speakers
    active = np.zeros((len(centers), len(speakers)), dtype=bool)
    for start, end, speaker, _ in clip.turns:
        active[:, speakers.index(speaker)] |= (centers >= start) & (centers < end)
    return active


def gss_waveforms(audio, clip):
    import numpy as np

    from codes.chapters.ch02.core.spectral import istft, stft
    from codes.chapters.ch07.core.dereverberation import offline_wpe
    from codes.chapters.ch08.core.gss_teaching import guided_cacgmm_mvdr

    n_fft, hop = 512, 160
    spectrum = stft(audio, n_fft=n_fft, hop_length=hop).transpose(1, 0, 2)
    dereverbed = offline_wpe(spectrum, taps=3, delay=2, iterations=2)
    active = activity(clip, np.arange(spectrum.shape[-1]) * hop)
    separated = [guided_cacgmm_mvdr(dereverbed, active, iterations=8, target_speaker=speaker)["output"]
                 for speaker in range(active.shape[1])]
    return istft(np.stack(separated), n_fft=n_fft, hop_length=hop, length=clip.samples)


def speaker_intervals(clip, speaker):
    intervals = []
    for start, end, label, _ in clip.turns:
        if label != speaker:
            continue
        if intervals and start <= intervals[-1][1]:
            intervals[-1] = (intervals[-1][0], max(intervals[-1][1], end))
        else:
            intervals.append((start, end))
    return intervals


def run_whisper(clips, directory, model_path, mode, *, use_history=True):
    import mlx.core as mx
    import mlx_whisper
    from mlx_whisper.transcribe import ModelHolder

    require_sha256(model_path / "weights.npz", MLX_WHISPER_SHA256)
    require_sha256(model_path / "config.json", MLX_WHISPER_CONFIG_SHA256)
    mx.set_default_device(mx.cpu)
    model = ModelHolder.get_model(str(model_path), mx.float32)
    mx.eval(model.parameters())

    def decode(waveform):
        return mlx_whisper.transcribe(waveform.astype("float32"), path_or_hf_repo=str(model_path),
                                      language="en", temperature=0.0, fp16=False,
                                      condition_on_previous_text=use_history)["text"]

    for clip in clips:
        channels = (1,) if mode == "mono-cuts" else CHANNELS
        audio = load_audio(clip, directory, channels)
        start = perf_counter()
        streams = [audio[0]] * len(clip.speakers) if mode == "mono-cuts" else gss_waveforms(audio, clip)
        if mode == "gss-full":
            decoded = [decode(stream) for stream in streams]
        else:
            decoded = [" ".join(decode(stream[first:last]) for first, last in speaker_intervals(clip, speaker))
                       for speaker, stream in zip(clip.speakers, streams, strict=True)]
        elapsed = perf_counter() - start
        yield {"id": clip.name, "streams": decoded, "seconds": elapsed,
               "pcm_sha256": {str(channel): PCM_SHA256[clip.name, channel] for channel in channels}}


@contextmanager
def dicow_model_class(model_path):
    # Make relative imports resolve to the eight verified upstream source files.
    with TemporaryDirectory(prefix="ami_dicow_") as scratch:
        directory = Path(scratch)
        for name, expected in DICOW_SOURCE_SHA256.items():
            code = (model_path / name).read_bytes()
            if sha256(code).hexdigest() != expected:
                raise ValueError(f"SHA-256 mismatch: {model_path / name}")
            (directory / name).write_bytes(code)
        package = ModuleType(directory.name)
        package.__path__ = [str(directory)]
        sys.modules[package.__name__] = package
        try:
            implementation = import_module(f"{package.__name__}.modeling_dicow")
            yield implementation.DiCoWForConditionalGeneration
        finally:
            for name in list(sys.modules):
                if name == package.__name__ or name.startswith(f"{package.__name__}."):
                    del sys.modules[name]


def decode_dicow_clip(model, processor, waveform, clip):
    import numpy as np
    import torch

    with torch.inference_mode():
        features = processor.feature_extractor(waveform, sampling_rate=SAMPLE_RATE, return_tensors="pt")
        frames = features.input_features.shape[-1] // 2
        a, b = activity(clip, np.arange(frames) * 320).T
        silence = ~(a | b)
        overlap = a & b
        # The published encoder expects silence, target solo, other solo, overlap.
        stno = np.stack(([silence, a & ~b, b & ~a, overlap], [silence, b & ~a, a & ~b, overlap]))
        input_features = features.input_features.repeat(2, 1, 1)
        attention_mask = torch.ones((2, frames * 2), dtype=torch.bool)
        ids = model.generate(input_features=input_features, attention_mask=attention_mask,
                             stno_mask=torch.as_tensor(stno, dtype=input_features.dtype),
                             language="en", task="transcribe", return_timestamps=True,
                             do_sample=False, num_beams=1)
        return processor.batch_decode(ids, skip_special_tokens=True)


def run_dicow(clips, directory, model_path):
    import torch
    from transformers import AutoProcessor

    require_sha256(model_path / "model.safetensors", DICOW_SHA256)
    torch.set_num_threads(4)
    processor = AutoProcessor.from_pretrained(model_path, local_files_only=True)
    with dicow_model_class(model_path) as implementation:
        model = implementation.from_pretrained(model_path, local_files_only=True, use_safetensors=True).eval()
        model.set_tokenizer(processor.tokenizer)
        for clip in clips:
            waveform = load_audio(clip, directory, (1,))[0]
            start = perf_counter()
            decoded = decode_dicow_clip(model, processor, waveform, clip)
            elapsed = perf_counter() - start
            yield {"id": clip.name, "streams": decoded, "seconds": elapsed,
                   "pcm_sha256": {"1": PCM_SHA256[clip.name, 1]}}


def score(clips, experiments):
    from meeteval.io import SegLST
    from meeteval.wer import cpwer
    from whisper.normalizers import EnglishTextNormalizer

    from codes.chapters.ch15.examples.score_meeting_transcripts import _summarize

    normalizer = EnglishTextNormalizer()
    ref = SegLST([
        {"session_id": clip.name, "speaker": speaker,
         "words": normalizer(" ".join(words for _, _, tag, words in clip.turns if tag == speaker))}
        for clip in clips for speaker in clip.speakers
    ])
    results = {}
    expected = {clip.name for clip in clips}
    lengths = {clip.name: clip.samples / SAMPLE_RATE for clip in clips}
    duration = sum(lengths.values())
    for name, path in experiments:
        if name in results:
            raise ValueError(f"Duplicate method name: {name}")
        rows = json.loads(Path(path).read_text(encoding="utf-8"))["clips"]
        ids = [row["id"] for row in rows]
        if len(ids) != len(expected) or set(ids) != expected:
            raise ValueError(f"{name}: need exactly one result for each reference excerpt")
        hypo = []
        times = {}
        for row in rows:
            seconds = row["seconds"]
            if isinstance(seconds, bool) or not isinstance(seconds, (float, int)) or not math.isfinite(seconds) or seconds <= 0:
                raise ValueError(f"{name}: invalid elapsed seconds for {row['id']}")
            streams = row["streams"]
            if not isinstance(streams, list) or any(not isinstance(text, str) for text in streams):
                raise ValueError(f"{name}: expected a list of speaker streams for {row['id']}")
            times[row["id"]] = seconds
            for speaker, words in enumerate(streams or [""]):
                hypo.append({"session_id": row["id"], "speaker": str(speaker), "words": normalizer(words)})
        by_clip = cpwer(ref, SegLST(hypo), reference_sort=False, hypothesis_sort=False)
        results[name] = {
            **_summarize(by_clip),
            "elapsed_seconds": sum(times.values()),
            "rtf": sum(times.values()) / duration,
            "by_clip": {key: {"errors": int(value.errors), "reference_words": int(value.length),
                              "seconds": times[key], "rtf": times[key] / lengths[key]}
                        for key, value in by_clip.items()},
        }
    return {"audio_seconds": duration, "results": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cutset", type=Path, required=True, help="Pinned ESPnet AMI development cutset")
    subparsers = parser.add_subparsers(dest="command", required=True)
    preparation = subparsers.add_parser("prepare", help="Make matching PCM excerpts from two AMI WAVs")
    preparation.add_argument("--original", type=Path, required=True)
    preparation.add_argument("--clips", type=Path, required=True)
    for name in ("gss-full", "gss-cuts", "mono-cuts", "dicow"):
        method = subparsers.add_parser(name, help=f"Run {name} from a local model snapshot")
        method.add_argument("--clips", type=Path, required=True)
        method.add_argument("--model", type=Path, required=True)
        if name != "dicow":
            method.add_argument("--no-history", action="store_true", help="Decode without previous-text conditioning")
    scoring = subparsers.add_parser("score", help="Score per-speaker text from any of the three runners")
    scoring.add_argument("--method", action="append", nargs=2, metavar=("NAME", "RESULT"), required=True)
    args = parser.parse_args()
    clips = load_clips(args.cutset)
    if args.command == "prepare":
        prepare(clips, args.original, args.clips)
        return
    if args.command == "score":
        output = score(clips, args.method)
    elif args.command == "dicow":
        output = {"clips": list(run_dicow(clips, args.clips, args.model))}
    else:
        output = {"clips": list(run_whisper(clips, args.clips, args.model, args.command,
                                           use_history=not args.no_history))}
    print(json.dumps(output, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
