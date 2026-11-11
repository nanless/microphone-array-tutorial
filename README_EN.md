# An Introduction to Microphone Arrays

**Hear where sound comes from.** Start with a recording and learn how a conference device finds a speaker, listens toward them, and handles echo.

**[Download the course introduction →](assets/readme/course-intro.html)** · [Chapter 1](chapters/01_problem-definition.md) · [Full contents](chapters/00_overview.md) · [PDF](dist/microphone-array-tutorial.pdf) · [简体中文](README.md)

<picture>
  <source media="(prefers-reduced-motion: reduce)" srcset="assets/readme/course-cover-still.svg">
  <img alt="Source A on the left makes a sound. The four microphones receive it in order from left to right. In the course introduction, move the source to the right to see the order reverse." src="assets/readme/course-cover.svg">
</picture>

## Start here

On the [course file page](assets/readme/course-intro.html), choose **Download raw file** and open the saved HTML in your browser. The Chinese course introduction includes the animation, three audio players, and the learning route. The demo and audio work offline without Python or a microphone.

Or listen here: **[One mic](assets/readme/lab/single.wav?raw=true) → [steered left](assets/readme/lab/left.wav?raw=true) → [steered right](assets/readme/lab/right.wav?raw=true)**. Turn down the volume first. Each clip is about three seconds. Long tones are on the left, pairs of short tones on the right, and all three clips use the same simulated recordings.

## Preface · Why record the same words four times?

Halfway through a call, someone at the far end of the table asks, “Can you hear me?” In your headphones you get half a sentence, some keystrokes, and an echo of what you just said. Those sounds pile up, but they give the device different clues.

Separate the four raw recordings. The voice reaches the nearest microphone first and the farthest one later. **To find the speaker**, slide the recordings and measure how far apart they are. **To listen toward them**, use that offset to line up the voice before averaging the channels. Your own words can echo back if the remote loudspeaker plays them and the microphones pick them up again. In that case the device has another clue: it knows what it just played.

This book starts with what you can hear. At each step, ask what signals you have, what they can tell you, and when they might lead you astray.

## Where to read next

| What you want to do | Where to start |
| --- | --- |
| Work out where a sound came from | [Ch. 1–4: The task, signals, geometry, and direction](chapters/01_problem-definition.md) |
| Listen toward someone or track them as they move | [Ch. 5: Beamforming](chapters/05_beamforming.md) · [Ch. 9: Tracking](chapters/09_source-tracking.md) |
| Handle an echo, a ringing room, or overlapping voices | [Ch. 6: Echo cancellation](chapters/06_aec.md) · [Ch. 7: Dereverberation](chapters/07_wpe-dereverberation.md) · [Ch. 8: Speech separation](chapters/08_speech-separation.md) |
| Put the methods on a device and test the results | [Ch. 10: Sync, calibration, and latency](chapters/10_engineering-practice.md) · [Ch. 11: Selection and validation](chapters/11_selection-guide.md) |

Start with Chapter 1 and consult the [math and notation appendix](chapters/14_appendix-symbols-math.md) when you need complex numbers or matrices. The [full contents](chapters/00_overview.md) also covers acoustic imaging, distributed arrays, and room simulation. The book is in Chinese.

**13 chapters · 2 appendices · 363 runnable exercises · 80 numbered figures**

[Code and audio by chapter](codes/chapters/README.md) · [Listening experiment explained](assets/readme/lab/README_EN.md) · [Implementation notes](codes/chapters/ch00/research/README.md) · [Web and PDF builds](scripts/README.md) · [Third-party sources and terms](codes/chapters/ch00/THIRD_PARTY.md)

The repository does not have a general root license. Third-party code and data retain their own terms.
