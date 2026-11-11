# An Introduction to Microphone Arrays

**One conference table, four microphones. How do you make the person at the far end come through clearly?**

[Start reading →](chapters/01_problem-definition.md) · [Learning route](#learning-route) · [PDF](dist/microphone-array-tutorial.pdf) · [简体中文](README.md)

<picture>
  <source media="(prefers-reduced-motion: reduce)" srcset="assets/readme/cover/room-still.svg">
  <img alt="A conference table from above. A person speaks on the left. Four microphones light up in order from nearest to farthest as the sound crosses the table." src="assets/readme/cover/room.svg">
</picture>

## Preface: four recordings of the same room

Halfway through a call, someone at the far end of the table asks, “Can you hear me?” In your headphones you hear half a sentence, a few keystrokes, and an echo of the last thing you said. Switch seats and the problem changes again.

Imagine listening to the four raw recordings separately. They all contain the same room, but the same words reach each microphone at a slightly different time. In the drawing, the speaker sits closest to mic 1. The direct sound reaches mics 2, 3 and 4 a little later. The differences can be hard to hear. Plot the four waveforms and you can see them.

![The same brief waveform shifts to the right from track 1 to track 4. Time runs to the right, and gold dots mark the starts. This is a schematic, not a measurement.](assets/readme/cover/four-tracks.svg)

You can use those delays in two ways. **To find the speaker**, measure how much a piece of sound shifts across recordings and use the microphone positions to infer a direction. That's **direction of arrival estimation**. **To listen toward the speaker**, use that direction's delays to align the voice before combining the recordings. That's **beamforming**. A keyboard in a different direction generally won't line up under those same delays. How much it weakens depends on frequency, microphone spacing, and the room.

The echo has a different clue. The far end's voice may play through the loudspeaker and feed back into the microphones. The device knows what it played, so it can use that signal to predict the echo. Room reflections, overlapping speakers, and a moving talker each pose their own problem. This book starts with the four recordings on the table and takes those problems apart one at a time.

**Listen now.** Two synthetic rhythms share the same recordings. [Hear one mic](assets/readme/lab/single.wav?raw=true), then [steer left](assets/readme/lab/left.wav?raw=true) and [steer right](assets/readme/lab/right.wav?raw=true). Turn down the volume first; each clip is about three seconds. [Try the interactive experiment →](assets/readme/lab/README_EN.md)

## Learning route

New to arrays? Start with Chapter 1. You can look at the figures and listen before you run code. Refer to the [math and notation appendix](chapters/14_appendix-symbols-math.md) when you need complex numbers, matrices, or decibels. If you have a specific problem, jump to the relevant chapter. The main text is in Chinese.

### 01　Find a direction in the recordings

[Ch. 1: What to keep](chapters/01_problem-definition.md) → [Ch. 2: How sound reaches the microphones](chapters/02_basics-signal-model.md) → [Ch. 3: How to place the mics](chapters/03_array-geometry.md) → [Ch. 4: Estimate direction from delays](chapters/04_doa-estimation.md)

Slide two waveforms to find the best match. Then work out whether the measured delay identifies one direction or several possible directions.

### 02　Point the pickup toward a speaker

[Ch. 5: Beamforming](chapters/05_beamforming.md) → [Ch. 9: Track a moving speaker](chapters/09_source-tracking.md)

Align and average the channels. Then steer the wrong way on purpose and listen to what happens to the target. If the speaker moves, use consecutive direction estimates together so that one bad estimate doesn't swing the result.

### 03　Tell echo, room tails, and overlapping voices apart

[Ch. 6: Cancel echo with the playback signal](chapters/06_aec.md) · [Ch. 7: Reduce room reverberation](chapters/07_wpe-dereverberation.md) · [Ch. 8: Separate overlapping speakers](chapters/08_speech-separation.md)

Start with what you hear: your own last sentence, a reflection trailing a word, or a second person talking over it? Each problem comes with different clues.

### 04　Test real recordings and put the methods on a device

[A stream recorded by sixteen mics](codes/chapters/ch02/real_audio/README.md) · [Ch. 10: Sync, calibration, and latency](chapters/10_engineering-practice.md) · [Ch. 11: System selection and validation](chapters/11_selection-guide.md)

Compare a single channel against the average in one real recording, and see what happens when the noise is correlated. Then set targets for your own device: working distance, acceptable latency, and how well the pickup should hold up in different conditions.

To draw a map of **where a machine is making noise**, see [Ch. 12: Acoustic imaging](chapters/12_acoustic-imaging.md). For microphones on separate devices with separate clocks, see [Ch. 13: Distributed arrays](chapters/13_distributed-enhancement.md).

## Book and companion material

[Full contents](chapters/00_overview.md) · [Code and audio by chapter](codes/chapters/README.md) · [Exercise and audio handbook](codes/chapters/ch00/research/05_exercises_and_audio.md) · [Room simulation](chapters/15_appendix-guide.md)

The book has **13 chapters, 2 appendices, 363 runnable exercises, and 80 numbered figures**, plus a figure for the room exercise. It distinguishes synthetic controls you can calculate by hand, room simulations, and real recordings so you can see what each result supports.

[Implementation research notes](codes/chapters/ch00/research/README.md) · [Third-party sources and terms](codes/chapters/ch00/THIRD_PARTY.md) · [Build instructions](scripts/README.md)

The repository does not have a general root license. Third-party code and data retain their own terms.
