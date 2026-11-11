# Microphone Array Signal Processing

**How does a group of microphones know whom to listen to?** Start with the order in which sound reaches four microphones. Work your way to a conference mic that listens toward a speaker and a robot that follows a moving voice.

[Hear the four mics in action](#listen-same-scene-two-directions) · [Beginner course route](#follow-the-course) · [Book PDF](dist/microphone-array-tutorial.pdf) · [简体中文](README.md)

<picture>
  <source media="(prefers-reduced-motion: reduce)" srcset="assets/readme/listening-room-still.svg">
  <img alt="Source A on the left reaches M1 first. Source B on the right reaches M4 first. The animation slows down the two arrival sequences." src="assets/readme/listening-room.svg">
</picture>

## Listen: same scene, two directions

Picture a table. Source A on the left makes a series of long tones. Source B on the right makes pairs of short tones. **Make a prediction: can you bring each side forward in turn without moving the microphones?**

1. **[▶ Single mic M1](assets/readme/lab/single.wav?raw=true)**: both the long tones and the pairs of short tones.
2. **[▶ Steer left to A](assets/readme/lab/left.wav?raw=true)**: the long tones stay; the pairs get quieter.
3. **[▶ Steer right to B](assets/readme/lab/right.wav?raw=true)**: this time the pairs come through.

Turn your device volume down first. The clips are mono, so a phone speaker works. If your browser downloads the WAV, open it in a player.

The changes come from applying different delays to the four channels. The clips share the same export gain. These synthetic sounds are chosen to make the contrast easy to hear.

**[Download the interactive lab on a computer](assets/readme/lab/array-lab.html)** → choose **Download raw file**, then open the saved HTML in a browser. No Python or mic permission needed. Steer toward A, then change B to the low tone. B stayed in the same place. Why can you hear it again?

The lab has an **English** switch.

## Preface: why make the first microphone wait?

Someone starts speaking at the far end of a conference table just as a keyboard starts clacking nearby. Average four recordings as they are and both sounds will usually still be there. What gives the array a chance is the **order of arrival**. Sounds on the left reach the left end first; sounds on the right reach the right end first. That timing can tell us where a sound came from and help us listen in that direction.

Say the same peak from A reaches the four mics at samples 0, 2, 4 and 6. Make M1 wait six samples, M2 wait four, and M3 wait two. M4 has already arrived last, so it needs no extra delay. All four peaks now land at sample 6.

![The horizontal axis counts samples. The hollow markers for M1, M2 and M3 move along the arrows to sample 6. M4 is already there. The four filled markers line up.](assets/readme/delay-and-sum.svg)

Now the average reinforces A: its peaks line up, so it keeps its original strength. B arrives in the opposite order. Apply the delays for A and B's peaks can meet its troughs, weakening the average. This is the simplest form of **beamforming**: change how the recordings are combined, and the array can favor a direction without moving.

Different directions alone do not guarantee suppression. Change B to **250 Hz** and steer toward A again: about **82% of B's single-mic power** remains. The low tone takes longer to complete a cycle, so the same timing differences shift it by only a small part of a cycle. That is the pattern for this book: find out why a method works, then find out when it falls short.

## Follow the course

Start at stop 1 or jump to the problem you have. The lab needs no math background. Later, [Appendix A](chapters/14_appendix-symbols-math.md) covers the trigonometry, complex numbers, and matrices used in the book. The chapters are written in Chinese.

### Stop 1 · Understand why the channels differ

[Ch. 1: The task](chapters/01_problem-definition.md) → [Ch. 2: Signal models](chapters/02_basics-signal-model.md) → [Ch. 3: Array geometry](chapters/03_array-geometry.md)

See which waveform arrives first and calculate how changing mic spacing affects the delay.

### Stop 2 · Locate a sound and listen toward it

[Ch. 4: Direction of arrival](chapters/04_doa-estimation.md) → [Ch. 5: Beamforming](chapters/05_beamforming.md)

Slide two recordings until they agree, then compare steering in the right and wrong directions.

### Stop 3 · Tell echo, reverberation, and overlapping speech apart

[Ch. 6: Echo cancellation](chapters/06_aec.md) · [Ch. 7: Dereverberation](chapters/07_wpe-dereverberation.md) · [Ch. 8: Speech separation](chapters/08_speech-separation.md)

Hearing your own voice on a call, hearing a room ring, and hearing two people talk at once are three different problems. Describe what you hear, then choose a method.

### Stop 4 · Put the methods on a device

[Ch. 9: Tracking](chapters/09_source-tracking.md) · [Ch. 10: Engineering](chapters/10_engineering-practice.md) · [Ch. 11: System selection](chapters/11_selection-guide.md)

Handle a bad direction estimate, check channel sync, and measure latency and pickup separately.

Prefer **a real recording**? Try [one stream recorded by sixteen microphones](codes/chapters/ch02/real_audio/README.md). Listen to the single mic and the average, then check why adding mics does not always cut noise power by the textbook ratio.

For acoustic maps or arrays spread across multiple devices, continue with [Ch. 12: Acoustic imaging](chapters/12_acoustic-imaging.md) and [Ch. 13: Distributed arrays](chapters/13_distributed-enhancement.md).

## Book and companion material

[Full contents](chapters/00_overview.md) · [Offline PDF](dist/microphone-array-tutorial.pdf) · [Code and audio by chapter](codes/chapters/README.md) · [Exercise and audio handbook](codes/chapters/ch00/research/05_exercises_and_audio.md)

The book has **13 chapters, 2 appendices, 363 runnable exercises, and 80 numbered figures**, plus a figure for the room simulation exercise. Use the explanations, hand calculations, code, and audio to check each other. [Appendix B](chapters/15_appendix-guide.md) includes a hands-on room simulation.

[Implementation research notes](codes/chapters/ch00/research/README.md) · [Third-party sources and terms](codes/chapters/ch00/THIRD_PARTY.md) · [Build instructions](scripts/README.md)

<details>
<summary>Reproduce the homepage experiment</summary>

The [lab HTML](assets/readme/lab/array-lab.html) contains the full synthesis and processing code. At 16 kHz the arrivals are `[0, 2, 4, 6]` samples for A and `[6, 4, 2, 0]` for B. Each percentage uses the known floating-point component for one source: its output sum of squares divided by its single-mic sum of squares, including the complete delay tail. The listening files sum both components and quantize the result to PCM WAV.

With Node.js installed, run `node assets/readme/lab/export-lab.mjs --check` from the repository root to regenerate the clips in memory and check the published bytes. The numbers describe this synthetic example; they are not measurements of a physical array.

</details>

The repository does not have a general root license. Third-party code and data retain their own terms.
