# 麦克风阵列入门

[English guide](README_EN.md) · [第一课：两支麦克风](chapters/01_problem-definition.md) · [全书导读](chapters/00_overview.md) · [下载 PDF](dist/microphone-array-tutorial.pdf)

**同一句话，两支麦克风都录到了，为什么一混合反而少了一部分声音？** 你在会议桌一端说话，近处的麦克风先收到，远处的晚一点。声源定位利用的正是这点先后，混合之前也得先想办法对齐它们。

![同一个声音沿不同长短的路径到达两支麦克风，右边的两条波形先后起振。](assets/readme/array-listening.svg)

## 先听：只差一个采样点，高音去哪儿了？

把一个低音和一个高音叠起来，让 M2 比 M1 晚记下一格。数字录音每隔固定时间记一个数，这一格就是一个**采样点**。下面是同一段合成声音的两种混合结果，没有加噪声，也没有分别调音量。先调低音量，按 A、B 各听一次，每段 2 秒。

**[▶ A · 直接把两路平均](codes/chapters/ch02/audio/alignment_unaligned.wav?raw=true)**　　**[▶ B · 先对齐，再平均](codes/chapters/ch02/audio/alignment_aligned.wav?raw=true)**

如果浏览器下载了 WAV，用本地播放器打开即可。比较时，留意 B 中较尖的那层音，低一点的音变化小得多。听不出也没关系：下一张图能告诉你变化来自哪里。

## 再看：同一次到达，为什么变成了两个小峰？

用一个短脉冲代替整段声音。图里的横轴是采样编号，每根竖线的高度就是这一点的振幅；`½ Σ` 表示把两路相加再除以二。

![脉冲对齐示意。A：两路分别在第 2、3 点出现峰值，直接平均留下两个半高峰。B：让先到的 M1 晚一格，两路在第 3 点叠出一个全高峰。](assets/readme/one-sample-shift.svg)

- **A，直接平均。** 第 2 点是 `(1 + 0) ÷ 2 = 0.5`，第 3 点是 `(0 + 1) ÷ 2 = 0.5`。原来的尖峰被分到了两格。
- **B，对齐再平均。** 让先到的 M1 等一格。两路的峰都落到第 3 点，得到 `(1 + 1) ÷ 2 = 1`。等的是先到那路，所以不需要预知下一点。

脉冲图解释了**怎么对齐**。刚才听的是连续的音，每一周期还有正、负两半。两路一旦错开，一个音的正值就可能碰上另一路的负值；平均后的振幅会变小。

## 算一步：同样晚一格，为什么高音损失更多？

这组音频每秒记录 16,000 点，一格是 **62.5 微秒**。1 kHz 的音每周期有 16 点，4 kHz 的音每周期只有 4 点。同样错开一格，前者只错开周期的 **1/16**，后者却错开了 **1/4**。

![从两段 WAV 读回的振幅。A 在 1 kHz 保留 98.1%，在 4 kHz 保留 70.7%；B 两个频率均保留 100%。](assets/readme/two-tone-result.svg)

图中的 100% 指[干净参照](codes/chapters/ch02/audio/alignment_reference.wav?raw=true)在相同频率上的振幅。4 kHz 直接平均后只剩约七成，1 kHz 则仍接近原值。这是**振幅**的比例；实验没有噪声，也没有测真人语音的清晰度。

<details>
<summary>想手算 70.7%？只用四个点</summary>

4 kHz 在 16 kHz 采样率下每周期正好四点。把原始正弦的振幅看作 1，一路依次是 `0, 1, 0, −1`；晚一格的那路是 `−1, 0, 1, 0`（取稳态的一周期）。

两路平均得到 `−0.5, 0.5, 0.5, −0.5`。这四个点碰巧都没踩在新正弦的峰上。要算峰值，可把平均写成 `0.5 × sin(θ) − 0.5 × cos(θ)`；这两个相差四分之一周期的分量合起来，振幅是 `√(0.5² + 0.5²) = 1/√2 ≈ 0.707`。完整解答和音频测量见 [E01-05](chapters/01_problem-definition.md#e01-05差一个采样点为什么高频损失更大)。

</details>

## 下一课，按你想解决的问题选

| 你想弄懂什么 | 从这里读 |
| --- | --- |
| **两支麦为什么有用？** | [第 1 章：哪些声音想保留，哪些想去掉](chapters/01_problem-definition.md) → [第 2 章：把到达先后写成模型](chapters/02_basics-signal-model.md) |
| **怎么知道声音从哪边来，又怎样朝那边拾音？** | [第 4 章：先找两路最吻合的错位](chapters/04_doa-estimation.md) → [第 5 章：让目标方向的声音对齐再合成](chapters/05_beamforming.md) |
| **为什么通话会听见自己的回声，空房间里又会拖尾？** | [第 6 章：借助设备知道的播放信号消回声](chapters/06_aec.md) · [第 7 章：处理房间反射造成的混响](chapters/07_wpe-dereverberation.md) |
| **要把算法装进设备，该测什么？** | [第 10 章：延迟、同步与标定](chapters/10_engineering-practice.md) · [第 11 章：从场景和验收条件反推选型](chapters/11_selection-guide.md) |

用到达先后推断方向，是**声源定位**的起点；按目标方向补时差、再合成，是最简单的**波束形成**。这里直接给了正确时差。真实录音里，如何找到它、什么时候会找错，正是后续课程要解决的问题。

<details>
<summary>在自己电脑上复算这段音频</summary>

需要 Git 和 Python。本书用 Python 3.13 验证；依赖版本见 [requirements.txt](requirements.txt)。

```bash
git clone https://github.com/nanless/microphone-array-tutorial.git
cd microphone-array-tutorial
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -c 'import json; from codes.chapters.ch01.chapter01_experiments import published_alignment_measurements as measure; print(json.dumps(measure()["pcm_measurements"], indent=2))'
```

Windows 将 `.venv/bin/python` 换成 `.venv\Scripts\python.exe`。结果分别列出 1 kHz 和 4 kHz 的振幅及比值。函数直接读取随仓库提供的 WAV，核对格式和摘要后再测量，不改写音频。

</details>

## 全书与资料

[章节导览](chapters/00_overview.md) · [离线合订 PDF](dist/microphone-array-tutorial.pdf) · [按章查代码与音频](codes/chapters/README.md) · [练习与音频手册](codes/chapters/ch00/research/05_exercises_and_audio.md)

全书含 **13 章正文、2 篇附录、363 道可运行练习和 80 张编号图**，另有一张房间题补图。可按问题跳读；几何、语音分离、追踪、声学成像和分布式阵列都在[全书导读](chapters/00_overview.md)中有各自入口。[音频目录](codes/chapters/README.md#audio-inventory)把合成实验与 DEMAND 真实录音分开列出。

[源码研究手册](codes/chapters/ch00/research/README.md) · [第三方来源与许可](codes/chapters/ch00/THIRD_PARTY.md) · [网站、PDF 和实验的构建说明](scripts/README.md)

仓库根目录没有通用许可证；公开可见不等于获准再分发。第三方代码和数据按各自条款使用。
