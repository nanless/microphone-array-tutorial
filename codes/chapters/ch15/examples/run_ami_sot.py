"""Decode the selected AMI excerpts with the ESPnet Whisper SOT release."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import time
import types

from codes.chapters.ch00.core.audio_samples import read_pcm16


MODEL_REVISION = "1b2452421d3acb2405955e66e8b668928632d9cb"
DECODER_REVISION = "9a2b4def1fce1c7b4be3b0d98483fd1d3142f928"
MODEL_SHA256 = "11008ed2cbf9399dadf3a31023c3a3e7e8daeb710cef4851dc6971a7522e3a07"
DECODER_SHA256 = "d758d2f0eaa46b1006f77a76fa89b2d1151032aa99b24d643ac676d50dc7b70e"
TIMESTAMP = re.compile(r"<\|[0-9]+\.[0-9]+\|>")


def _load_decoder(path):
    code = path.read_bytes()
    if hashlib.sha256(code).hexdigest() != DECODER_SHA256:
        raise ValueError("expected the fixed ESPnet AMI decoder source")
    decoder = types.ModuleType("espnet_ami_sot_decoder")
    exec(compile(code, str(path), "exec"), vars(decoder))
    return decoder


def _streams(result, tokenizer, separator):
    text = "".join(tokenizer.decode_with_timestamps(segment["tokens"])
                   for segment in result["segments"])
    return [clean for block in text.split(separator)
            if (clean := " ".join(TIMESTAMP.sub("", block).split()))]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--decoder", required=True, type=Path, help="official local/decode.py")
    parser.add_argument("--checkpoint", required=True, type=Path, help="official model.pth")
    parser.add_argument("clips", nargs="+", type=Path, help="16 kHz mono PCM16 WAV files")
    args = parser.parse_args()

    with args.checkpoint.open("rb") as stream:
        if hashlib.file_digest(stream, "sha256").hexdigest() != MODEL_SHA256:
            raise ValueError("expected the fixed ESPnet AMI model checkpoint")
    decoder = _load_decoder(args.decoder)
    import torch
    import whisper

    torch.set_num_threads(4)
    tokenizer = whisper.tokenizer.get_tokenizer(multilingual=True, task="transcribe", language="en")
    separator = decoder.DEFAULT_SEPARATOR
    separator_tokens = tokenizer.encode(separator)
    if len(separator_tokens) != 1:
        raise ValueError("the ESPnet separator must be a single Whisper token")
    decoder._install_sot_separator_patch(separator_tokens[0])

    # The checkpoint supplies every parameter. Only the base architecture is needed.
    dimensions = whisper.model.ModelDimensions(
        n_mels=80, n_audio_ctx=1500, n_audio_state=768, n_audio_head=12, n_audio_layer=12,
        n_vocab=51865, n_text_ctx=448, n_text_state=768, n_text_head=12, n_text_layer=12,
    )
    state = torch.load(args.checkpoint, map_location="cpu", weights_only=True)
    model = whisper.model.Whisper(dimensions).eval()
    model.load_state_dict(decoder.remap_espnet_to_whisper(state), strict=True)
    del state

    clips = []
    for path in args.clips:
        wav = path.read_bytes()
        rate, audio = read_pcm16(wav)
        if rate != 16_000 or audio.shape[0] != 1 or not 0 < audio.shape[1] <= 30 * rate:
            raise ValueError(f"expected a mono 16 kHz WAV up to 30 seconds: {path.name}")
        waveform = audio[0].astype("float32")
        start = time.perf_counter()
        result = model.transcribe(
            waveform, language="en", task="transcribe", beam_size=5, fp16=False,
            word_timestamps=False, condition_on_previous_text=False, max_initial_timestamp=None,
        )
        clips.append({"id": path.stem.split(".")[0], "wav_sha256": hashlib.sha256(wav).hexdigest(),
                      "streams": _streams(result, tokenizer, separator),
                      "audio_seconds": len(waveform) / rate, "seconds": time.perf_counter() - start})

    print(json.dumps({"model_revision": MODEL_REVISION, "model_sha256": MODEL_SHA256,
                      "decoder_revision": DECODER_REVISION, "decoder_sha256": DECODER_SHA256,
                      "torch": torch.__version__, "whisper": whisper.__version__,
                      "device": "cpu", "torch_threads": torch.get_num_threads(), "clips": clips},
                     indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
