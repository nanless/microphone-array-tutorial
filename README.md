# 麦克风阵列入门

**把声音听出方向。** 从一段录音学起，理解会议机怎样找人、朝人拾音，以及为什么有些回声得另想办法。

**[在线进入课程 →](https://nanless.github.io/microphone-array-tutorial/)** · [直接读第一章](chapters/01_problem-definition.md) · [全书目录](chapters/00_overview.md) · [PDF](dist/microphone-array-tutorial.pdf) · [English](README_EN.md)

<picture>
  <source media="(prefers-reduced-motion: reduce)" srcset="assets/readme/course-cover-still.svg">
  <img alt="左边 A 发声，四支麦克风由左到右依次收到声音。课程页可以把发声位置换到右边，观察次序如何反转。" src="assets/readme/course-cover.svg">
</picture>

## 第一次来，从这里开始

[打开在线课程](https://nanless.github.io/microphone-array-tutorial/)。点选左右声源看动图，用页面上的三个播放器听差别，再按学习路线选章节。手机和电脑都能直接打开。

也可以就在 GitHub 上先听：**[单支麦](assets/readme/lab/single.wav?raw=true) → [对准左边](assets/readme/lab/left.wav?raw=true) → [对准右边](assets/readme/lab/right.wav?raw=true)**。先调低音量。每段约三秒，左边是一串长音，右边是成对短音，三段用的是同一组模拟录音。

## 序言 · 同一句话，为什么录了四份？

会开到一半，桌子那头问：「你们听得见吗？」耳机里传来半句话、键盘声，还有自己刚说过的话。听起来都像「声音太杂」，却是三道不同的题。

先把会议机里的四路录音分开。说话声先到近处的麦克风，晚一点才到远处。**找人时**，把两条波形前后滑动，看它们错开了多少；**拾音时**，反过来用这个时差把那个人的话对齐，再取平均。自己说过的话之所以绕回来，可能是它从对方的扬声器播出、又被麦克风收了进去。处理它还有另一条线索：设备知道自己刚播了什么。

这本书从这些听得见的现象出发。每一步都问：手上有什么信号，能靠它知道什么，什么时候会猜错？

## 这本书怎么读

| 手头的问题 | 从这里进去 |
| --- | --- |
| 第一次学，想从两路录音看出方向 | [1–4 章：任务、信号、几何、定位](chapters/01_problem-definition.md) |
| 找到人了，想让设备朝他拾音 | [第 5 章：波束形成](chapters/05_beamforming.md) · [第 9 章：跟踪移动声源](chapters/09_source-tracking.md) |
| 听见回声、房间拖尾，或两个人同时说话 | [第 6 章：消回声](chapters/06_aec.md) · [第 7 章：去混响](chapters/07_wpe-dereverberation.md) · [第 8 章：语音分离](chapters/08_speech-separation.md) |
| 装到设备上，想知道怎么测才算有效 | [第 10 章：同步、标定和延迟](chapters/10_engineering-practice.md) · [第 11 章：选型与验收](chapters/11_selection-guide.md) |

从第一章开始读，不用先学完复数和矩阵；遇到时再翻[数学与符号附录](chapters/14_appendix-symbols-math.md)。[全书目录](chapters/00_overview.md)另有声学成像、分布式阵列和房间仿真路线。

**13 章正文 · 2 篇附录 · 363 道可运行练习 · 80 张编号图**

[按章找代码和音频](codes/chapters/README.md) · [试听实验的推导](assets/readme/lab/README.md) · [研究手册](codes/chapters/ch00/research/README.md) · [网页与 PDF 构建](scripts/README.md) · [第三方来源与许可](codes/chapters/ch00/THIRD_PARTY.md)

仓库根目录没有通用许可证。第三方代码和数据按各自条款使用。
