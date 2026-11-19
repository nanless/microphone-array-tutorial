# 麦克风阵列信号处理 · 从入门到工业实践

从三角几何、复相量时频、矩阵优化与概率统计预备出发，结合多路录音和阵列几何，逐步推导声源定位、波束形成、回声消除、去混响、分离和追踪，再进入设备实现、评测、选型与验收。书中另有声学成像、分布式阵列、空间可穿戴阵列与条件目标提取四篇专题。

**[在线进入课程 →](https://nanless.github.io/microphone-array-tutorial/)** · [数学预备（第 0 章）](chapters/00_math-foundations.md) · [从第一章开始](chapters/01_problem-definition.md) · [工业实践](https://nanless.github.io/microphone-array-tutorial/site/#course-engineering) · [全书目录](chapters/00_overview.md) · [PDF](dist/microphone-array-tutorial.pdf) · [English](README_EN.md)

<picture>
  <source media="(prefers-reduced-motion: reduce)" srcset="assets/readme/course-cover-still.svg">
  <img alt="左边 A 发声，四支麦克风由左到右依次收到声音。课程页可以把发声位置换到右边，观察次序如何反转。" src="assets/readme/course-cover.svg">
</picture>

## 在线案例

[打开在线课程](https://nanless.github.io/microphone-array-tutorial/)。可用动图看声音到达的先后，再用同一段房间实录看定位、听原始录音、对齐相加和 MVDR 加后滤波三种输出。 [直接看定位与波束形成案例 →](https://nanless.github.io/microphone-array-tutorial/site/#course-cases)

## 序言 · 同一句话，为什么录了四份？

会开到一半，桌子那头问：「你们听得见吗？」耳机里传来半句话、键盘声，还有自己刚说过的话。听起来都像「声音太杂」，却是三道不同的题。

先把会议机里的四路录音分开。说话声先到近处的麦克风，晚一点才到远处。**找人时**，把两条波形前后滑动，看它们错开了多少。

**拾音时**，反过来用这个时差把那个人的话对齐，再取平均。自己说过的话之所以绕回来，可能是它从对方的扬声器播出、又被麦克风收了进去。处理它还有另一条线索：设备知道自己刚播了什么。

这本书从这些听得见的现象出发。每一步都问：手上有什么信号，能靠它知道什么，什么时候会猜错？

## 学习路线

| 阶段 | 可以做什么 | 从哪读起 |
| --- | --- | --- |
| 第 0～3 章 · 数学预备与阵列几何 | 夯实复相量、子空间投影、Wirtinger 微积分与贝叶斯工具，用时差、相位和阵列摆位建立信号模型 | [第 0 章 · 数学基础](chapters/00_math-foundations.md) · [第 1 章](chapters/01_problem-definition.md) |
| 第 4～9 章 · 算法设计与验证 | 按任务选择定位、波束形成、消回声、去混响、分离或追踪，复算并运行例子 | [声源定位](chapters/04_doa-estimation.md) · [波束形成](chapters/05_beamforming.md) · [算法目录](chapters/00_overview.md) |
| 第 10～11 章 · 工业实践 | 做同步和标定，测端到端延迟，核算算力与功耗，制定选型和验收指标 | [工程实现](chapters/10_engineering-practice.md) · [选型与验收](chapters/11_selection-guide.md) |
| 第 12～15 章 · 专题 | 定位设备噪声源、跨设备窄带协同、可穿戴空间拾音，或做条件目标提取与多通道日志 | [声学成像](chapters/12_acoustic-imaging.md) · [分布式阵列](chapters/13_distributed-enhancement.md) · [可穿戴阵列](chapters/14_spatial-wearable-arrays.md) · [目标提取与日志](chapters/15_target-extraction-diarization.md) |

第一次学可以先通读或按需检索第 0 章数学预备，也可以直接从第 1 章建立物理直觉，遇到具体符号与紧凑公式再速查附录 A。正在做设备，可以直接进入工程任务和工业源码手册。 [第 0 章数学基础](chapters/00_math-foundations.md) · [符号与公式速查](chapters/16_appendix-symbols-math.md) · [工程任务入口](https://nanless.github.io/microphone-array-tutorial/site/#course-engineering) · [工业源码手册](codes/chapters/ch00/research/03_industrial_deployment.md)

**16 章正文 · 2 篇附录 · 443 道可运行练习 · 92 张编号图**

[按章找代码和音频](codes/chapters/README.md) · [研究手册](codes/chapters/ch00/research/README.md) · [网页与 PDF 构建](scripts/README.md) · [第三方来源与许可](codes/chapters/ch00/THIRD_PARTY.md)

仓库根目录没有通用许可证。第三方代码和数据按各自条款使用。
