"""历史标题锚：章节重排后仍指向原来的教学主题。

键是旧版发布过的片段，值是当前标题片段。改版前的序号别名由固定
清单保存；语义编号和无编号标题摘要的少量变化也在本文件登记。
"""

import json
from pathlib import Path

try:
    from scripts.chapter_identity import historical_source, numbering_heading_aliases
except ModuleNotFoundError:
    from chapter_identity import historical_source, numbering_heading_aliases


# 改版前发布过的 sec-N 指向当时第 N 个标题。新标题增删后不能把
# 它机械地重新指向当前第 N 个标题，故固定旧版逐标题映射。
HISTORICAL_SEQUENTIAL_IDS = json.loads(
    Path(__file__).with_name("legacy_sequential_anchors.json").read_text(encoding="utf-8")
)

HISTORICAL_SECTION_IDS = {
    "14_acoustic-imaging.md": {
        "sec-u-c93fb52739": "sec-u-8378eb57b1",
        "sec-u-d6bcd65a94": "sec-u-8c3fd50486",  # 三次前向sweep -> 三次前向遍历
        "sec-u-cdd0dc66e7": "sec-u-9936ab612e",  # E12-06：控制 -> 基准对照
        "sec-u-f0f86235ee": "sec-12-13",         # old E12-02 pruned into prose
        "sec-u-1844251c9d": "sec-u-6afcc47205",  # old E12-03 -> E12-02
        "sec-u-ab20495b47": "sec-u-10593ce475",  # old E12-04 -> E12-03
        "sec-u-2fcef4767e": "sec-u-5f859559c9",  # old E12-05 -> E12-04
        "sec-u-9936ab612e": "sec-u-0295e9517e",  # old E12-06 -> E12-05
        "sec-u-e1197de104": "sec-u-fc79dc44f2",  # old E12-07 -> E12-06
        "sec-u-864e5d6a6f": "sec-12-13",         # old E12-08 pruned
        "sec-u-7731a0181a": "sec-u-25be8b6461",  # old E12-09 -> E12-07
        "sec-u-c1484bb4f1": "sec-u-b9276fe765",  # old E12-10 -> E12-08
        "sec-u-ce6e25be4f": "sec-u-95831eb566",  # old E12-11 -> E12-09
        "sec-u-f72f530d8f": "sec-u-fdc10bae1e",  # old E12-12 -> E12-10
        "sec-u-54ce6bf1ed": "sec-u-9b2b724fdf",  # old E12-13 -> E12-11
        "sec-u-c010f17b9c": "sec-u-df2535e30e",  # old E12-14 -> E12-12
        "sec-u-ca5815de33": "sec-u-65dbe997d4",  # old E12-15 -> E12-13
        "sec-u-3a315b135c": "sec-u-6454fae4c0",  # old E12-16 -> E12-14
        "sec-u-35150a13e5": "sec-u-cd2132f69a",  # old E12-17 -> E12-15
        "sec-u-21b3677ebd": "sec-u-25728c151b",  # old E12-18 -> E12-16
    },
    "15_distributed-enhancement.md": {
        "sec-u-d261d70382": "sec-u-89268d58ca",  # 作者MATLAB静态合同 -> 工程实现细节
        "sec-u-89268d58ca": "sec-u-7b11ea0af8",  # 工程实现细节 -> 在线流式实现的平滑延迟、窗归一化与门控对齐
        "sec-u-b345641c53": "sec-u-a6af0cfe91",  # 可听控制 -> 可听对照音频
        "sec-u-320ff7d3af": "sec-u-87bec6abd7",  # 采样率偏差控制 -> 采样率偏差补偿
        "sec-u-ae60fe614c": "sec-u-c8f4d1ccc4",  # E13-01 全角空格 -> 冒号
        "sec-u-d2277f09c0": "sec-u-b01a879a4a",  # E13-02
        "sec-u-990fd59df5": "sec-u-db8a2119c9",  # E13-03
        "sec-u-9493aad880": "sec-u-ada0f6befb",  # E13-04
        "sec-u-a123bcf090": "sec-u-275e8f69e6",  # E13-05
        "sec-u-b7dac779da": "sec-u-0f09bd925b",  # E13-06
        "sec-u-c8f4d1ccc4": "sec-13-16",         # old E13-01..06 folded into §13.16
        "sec-u-b01a879a4a": "sec-13-16",
        "sec-u-db8a2119c9": "sec-13-16",
        "sec-u-ada0f6befb": "sec-13-16",
        "sec-u-275e8f69e6": "sec-13-16",
        "sec-u-0f09bd925b": "sec-13-16",
        "sec-u-9bcc32b7fe": "sec-u-66bc517dc7",  # old E13-07 -> E13-01
        "sec-u-66bc517dc7": "sec-u-e712089bd9",
        "sec-u-40ecac23db": "sec-u-86e56e7136",  # old E13-08 -> E13-02
        "sec-u-86e56e7136": "sec-u-1339e48bd7",
        "sec-u-15a9194ba1": "sec-u-2a988d60f2",  # old E13-09 -> E13-03
        "sec-u-2a988d60f2": "sec-u-745773ac58",
        "sec-u-ecaa4b04e4": "sec-u-76c29b9e9a",  # old E13-10 -> E13-04
        "sec-u-76c29b9e9a": "sec-u-243b80151f",
        "sec-u-5159378837": "sec-u-8f41df6077",  # old E13-11 -> E13-05
        "sec-u-8f41df6077": "sec-u-16ab690f45",
        "sec-u-481fa3de3e": "sec-u-1ba11b3526",  # old E13-12 -> E13-06
        "sec-u-1ba11b3526": "sec-u-ee42fea78d",
        "sec-u-43c36175b8": "sec-u-55891f9c74",  # old E13-13 -> E13-07
        "sec-u-55891f9c74": "sec-u-01d67224fd",
        "sec-u-478a33071c": "sec-u-f54078f6df",  # old E13-14 -> E13-08
        "sec-u-f54078f6df": "sec-u-75575e4f83",
        "sec-u-4bad0eb87f": "sec-u-7cead60afd",  # old E13-15 -> E13-09
        "sec-u-7cead60afd": "sec-u-7ec15b4413",
        "sec-u-d6bfdb5b5e": "sec-u-9f6977a658",  # old E13-16 -> E13-10
        "sec-u-9f6977a658": "sec-u-7af8af3cc0",
        "sec-u-270f9db859": "sec-u-2ab25b15a4",  # old E13-17 -> E13-11
        "sec-u-2ab25b15a4": "sec-u-ce0e0d840a",
        "sec-u-4e97b65795": "sec-u-03e1ca65d2",  # old E13-18 -> E13-12
        "sec-u-03e1ca65d2": "sec-u-97a3c2617d",
        "sec-u-0c9a031714": "sec-u-87f37a1d30",  # old E13-19 -> E13-13
        "sec-u-87f37a1d30": "sec-u-826f17e30f",
        "sec-u-5bc0cc2c13": "sec-u-0dcd9a42fe",  # old E13-20 merged into E13-13
        "sec-u-0dcd9a42fe": "sec-u-826f17e30f",
        "sec-u-9e15aaffec": "sec-u-2b170d53fc",  # old E13-21 -> E13-14
        "sec-u-2b170d53fc": "sec-u-641c7c2c6a",
        "sec-u-4db262b2cd": "sec-u-0469bde73b",  # old E13-22 -> E13-15
        "sec-u-0469bde73b": "sec-u-3631413245",
        "sec-u-0b8e41071b": "sec-u-c7f8fcc136",  # old E13-23 -> E13-16
        "sec-u-c7f8fcc136": "sec-u-d07c497049",
        "sec-u-7e45f79bb0": "sec-u-f7034ff89c",  # old E13-24 -> E13-17
        "sec-u-f7034ff89c": "sec-u-63cbc678eb",
        "sec-u-632be6e3e5": "sec-u-13a2e4109a",  # old E13-25 -> E13-18
        "sec-u-13a2e4109a": "sec-u-c4927d5b7e",
        "sec-u-60abc8d82d": "sec-u-1979c49802",  # old E13-26 -> E13-19
        "sec-u-1979c49802": "sec-u-77723688f3",
        "sec-u-1f2fb8b5ba": "sec-u-8f5d54d722",  # old E13-27 -> E13-20
        "sec-u-8f5d54d722": "sec-u-45389cc58a",
    },
    "16_spatial-wearable-arrays.md": {
        "sec-u-82e576bd7f": "sec-u-bfe748c5ee",  # 锚点手算 14-3 图 82 -> 图 77
        "sec-u-2446522ba9": "sec-u-f3056f87a9",  # 锚点手算 14-4（下）图 84 -> 图 79
        "sec-u-4d29765748": "sec-u-9f4b3966bf",  # old E14-05 -> 图 76 解读
        "sec-u-3bbb9fdf58": "sec-u-6b1501756c",  # old E14-06 -> E14-05
        "sec-u-998d17d4b8": "sec-u-c5f9a35654",  # old E14-07 -> E14-06
        "sec-u-7932ba5f44": "sec-u-c5f9a35654",  # old E14-08 -> E14-06
        "sec-u-2a0c89082a": "sec-u-b6d6b23c09",  # old E14-09 -> E14-07
        "sec-u-1772a2308a": "sec-u-e70b9483ab",  # old E14-10 -> 双耳 MVDR 共线坍缩
        "sec-u-e3dbfefc58": "sec-u-bfb7be7dca",  # old E14-11 -> E14-08
        "sec-u-4b9298cdaa": "sec-u-e47db953e2",  # old E14-12 -> E14-09
        "sec-u-6385bcdab7": "sec-u-b498770aad",  # old E14-13 -> E14-10
        "sec-u-6c58daa8de": "sec-u-40f9972c5d",  # old E14-14 -> E14-11
        "sec-u-d9fe9adc55": "sec-u-8463384ecc",  # old E14-15 -> E14-12
        "sec-u-5a34d829b8": "sec-u-4ed81e38ff",  # old E14-16 -> E14-13
        "sec-u-a2a77f50a2": "sec-u-52f6b499cb",  # old E14-17 -> E14-14
        "sec-u-d3076b1ca1": "sec-u-66209d9215",  # old E14-18 -> E14-15
        "sec-u-1e132e185f": "sec-u-318398f81f",  # 扩散度柯西—施瓦茨 -> 严格界证明
    },
    "17_target-extraction-diarization.md": {
        "sec-u-7260e8bb52": "sec-15-10",         # old E15-01 pruned
        "sec-u-477f784bc2": "sec-u-ad12f37af4",  # old E15-02 -> E15-01
        "sec-u-972638914c": "sec-u-2d3d839fde",  # old E15-03 -> E15-02
        "sec-u-dba20ea272": "sec-15-3",          # old E15-04 -> §15.3 近场球面波
        "sec-u-3219c8f9e8": "sec-u-c8721eee93",  # old E15-05 -> E15-03
        "sec-u-f8d6da3847": "sec-u-3f00ade651",  # old E15-06 -> E15-04
        "sec-u-c101fdacd3": "sec-u-81ff5862ba",  # old E15-07 -> E15-05
        "sec-u-cba35fb5e1": "sec-u-7fc253decf",  # old E15-08 -> E15-06
        "sec-u-babc01c78f": "sec-u-7492136e0e",  # old E15-09 -> E15-07
        "sec-u-98c35aede6": "sec-u-555b43ef15",  # old E15-10 -> E15-08
        "sec-u-f030237235": "sec-u-a724bb1ad5",  # old E15-11 -> E15-09
        "sec-u-76e2b3cbad": "sec-u-5cd392eb9e",  # old E15-12 -> E15-10
        "sec-u-5427085a40": "sec-15-6",          # old E15-13 -> §15.6 帧率异构
        "sec-u-fcb2ba8992": "sec-u-9e65f470a2",  # old E15-14 -> E15-11
        "sec-u-617b8246f6": "sec-15-6",          # old E15-15 -> §15.6 唇部遮挡
        "sec-u-22086f7ea6": "sec-u-3e004b09fa",  # old E15-16 -> E15-12
        "sec-u-c7cc19996d": "sec-u-a6a7f14956",  # old E15-17 -> E15-13
        "sec-u-6fb2bbba75": "sec-u-a423afa58c",  # old E15-18 -> E15-14
    },
    "00_overview.md": {
        "sec-u-45eb1ba17e": "sec-u-bf0d4d6877",  # 92→87张图，同一插图地图
        "sec-u-d200046a3a": "sec-u-45eb1ba17e",  # 88→92张图，同一插图地图
        "sec-u-142ea3f046": "sec-u-a8458c949f",  # 18文件→19文件，同一文件地图
        "sec-u-bd5fecaa8e": "sec-u-d200046a3a",  # 80→88张图，同一插图地图
        "sec-u-461c9549d0": "sec-u-142ea3f046",  # 16文件→18文件，同一文件地图
        "sec-u-d12b48e159": "sec-u-bd5fecaa8e",  # 79→80张图，同一插图地图
        "sec-u-10c3704f3f": "sec-u-d12b48e159",  # 78→79张图，保留原插图地图主题
        "sec-u-f96924fd43": "sec-u-10c3704f3f",  # 77→78张图，同一插图地图
        "sec-u-8838edd879": "sec-u-f96924fd43",  # 76→77张图，仍是原插图地图
        "sec-u-55e537d650": "sec-u-8838edd879",  # 75→76张图，保留同一主题链
        "sec-u-dccee38a3d": "sec-u-55e537d650",  # 74→75张图，仍是原插图地图
        "sec-u-e445ebf98d": "sec-u-dccee38a3d",  # 73→74张图，仍是原插图地图
        "sec-u-726cf32300": "sec-u-e445ebf98d",  # 72→73张图，仍是同一插图地图
        "sec-u-dacd520d2c": "sec-u-726cf32300",  # 明示编号图，仍是同一插图地图
        "sec-u-2afd7c6509": "sec-u-461c9549d0",  # 15文件→16文件，同一文件地图
        "sec-u-9c77d218fe": "sec-u-dacd520d2c",  # 69图→72图，同一插图地图
        "sec-u-06e6cf7b78": "sec-u-2afd7c6509",  # 14文件→15文件，原文件地图
        "sec-u-ec3f837100": "sec-u-9c77d218fe",  # 66图→69图，仍是同一插图地图
        "sec-u-8657cda209": "sec-u-ec3f837100",  # 65 figure map keeps its original topic
        "sec-u-0d7a1562fa": "sec-u-ec3f837100",  # 64 张历史标题
        "sec-u-8a5dd76e40": "sec-u-ec3f837100",  # 63 张历史标题
        "sec-u-f19f47c7c1": "sec-u-ec3f837100",  # 61 张历史标题
        "sec-u-2b17484f3b": "sec-u-ec3f837100",  # 60 张历史标题
        "sec-u-cb1ede6c60": "sec-u-ec3f837100",  # 59 张历史标题
        "sec-u-5edae53c5f": "sec-u-ec3f837100",  # 58 张历史标题
        "sec-u-5c63f4d83e": "sec-u-ec3f837100",  # 56 张历史标题
        "sec-u-36fa20efde": "sec-u-ec3f837100",  # 54 张历史标题
        "sec-u-3ebc13ca5a": "sec-u-ec3f837100",  # 52 张历史标题
        "sec-u-1e72416790": "sec-u-ec3f837100",  # 50 张历史标题
        "sec-u-3a0278b879": "sec-u-ec3f837100",  # 49 张历史标题
        "sec-u-b7a71b077a": "sec-u-3937b1b94e",  # 路径 A：去掉无依据的周数
        "sec-u-efc5552983": "sec-u-f958564d39",  # 路径 B：去掉无依据的周数
        "sec-u-6ef8e18126": "sec-u-ec3f837100",  # 插图地图：40 张历史标题
        "sec-u-d538d6d0a5": "sec-u-ec3f837100",  # 41 张历史标题
        "sec-u-0943a9ed3c": "sec-u-ec3f837100",  # 43 张历史标题
        "sec-u-1e5d8a2bad": "sec-u-ec3f837100",  # 44 张历史标题
        "sec-u-d8fd1002de": "sec-u-ec3f837100",  # 45 张历史标题
        "sec-u-627ed3907e": "sec-u-ec3f837100",  # 47 张历史标题
    },
    "01_problem-definition.md": {
        "sec-u-81896116fa": "sec-u-5e074a7e92",  # 理想收益标题，原任务与收益主题
        "sec-u-5e074a7e92": "sec-1-1",           # 1.1 任务定义与理想阵列收益
        "sec-u-05ec932a34": "sec-u-c36b7ac16e",  # 客厅里的远场语音问题：标题改为远距离
        "sec-u-c36b7ac16e": "sec-1-1-1",         # 1.1.1 客厅里的远距离语音问题
        "sec-u-1d5f9c4fa9": "sec-1-2",           # 1.2 双耳听觉的三类定位线索
        "sec-u-f7a201edb9": "sec-1-3",           # 1.3 本章练习
    },
    "04_doa-estimation.md": {
        "sec-u-442f837dd2": "sec-u-5d3aebe0fb",  # MUSIC 旧标题
        "sec-u-5d3aebe0fb": "sec-u-cb730f7339",  # MUSIC、root-MUSIC 与前后向空间平滑
        "sec-u-8dda6571c9": "sec-u-e9bdae6530",  # ESPRIT 旧标题
        "sec-u-e9bdae6530": "sec-u-fd14b8a21b",  # ESPRIT（LS、TLS 与圆阵相位模式）
        "sec-u-1e0973b8ca": "sec-u-eb21d54fe3",  # 算例 4-3 R -> R_{xx}
        "sec-4-1-1": "sec-4-6-1",                # 原 4.1.1 源数估计 -> 4.6.1
        "sec-u-08146bdaaf": "sec-u-6fc8fe5b4a",  # 宽带聚焦：CSSM、RSS、WAVES 与 TOPS
        "sec-4-7-1": "sec-u-00305faff5",         # 近场模型失配
        "sec-4-7-2": "sec-u-374bae76e9",         # 可执行基线与外部实现
        "sec-u-374bae76e9": "sec-u-5e48a77123",  # 可执行基线与工程实现注意事项
        "sec-u-467e2bd676": "sec-4-8",           # 重复抽样 -> 4.8 定位算法总对比与分辨率边界
        "sec-4-10": "sec-u-58a3f5a5c9",          # 原 4.10 四麦亚采样时差并入 §4.2
        "sec-4-11": "sec-4-9",                   # 原 4.11 本章练习 -> 4.9
    },
    "05_beamforming.md": {
        "sec-u-fb93744989": "sec-5-5-1",  # 多通道维纳滤波与秩一目标模型 -> 5.5.1
        "sec-u-b11411e989": "sec-5-5-2",  # Frost 波束形成器：时域约束更新 -> 5.5.2
    },
    "03_array-geometry.md": {
        "sec-u-5834a7d67e": "sec-3-3-1",
        "sec-u-20a0d88f74": "sec-3-3-2",
        "sec-u-d14b5899a5": "sec-3-3-3",
        "sec-u-1d1a0d7440": "sec-3-1-2",  # 算例 3-1 -> 3.1.2
        "sec-u-83d5019885": "sec-3-1-4",  # 算例 3-2 -> 3.1.4
    },
    "06_aec.md": {
        "sec-6-1-3": "sec-6-2",
        **{f"sec-6-1-{index}": f"sec-6-{index - 1}"
           for index in range(4, 19)},
        "sec-6-11": "sec-6-8-2",  # 测试覆盖与验收矩阵 -> 6.8.2
        "sec-6-12": "sec-6-8-3",  # AEC 诊断方法与 NLP -> 6.8.3
        "sec-6-13": "sec-6-5-3",  # 个性化 AEC 合规 -> 6.5.3
        "sec-6-14": "sec-6-5-2",  # 研究方向 -> 6.5.2
        "sec-6-15": "sec-6-8-4",  # 可执行基线与开源边界 -> 6.8.4
        "sec-6-16": "sec-6-9",    # 常见误区与自测题 -> 6.9
        "sec-6-17": "sec-6-10",   # 延伸阅读 -> 6.10
        "sec-u-15fd0add02": "sec-u-920aee2313",  # 参考取点前移 -> 后移
    },
    "07_wpe-dereverberation.md": {
        "sec-7-1-1": "sec-7-8",
        "sec-7-1-2": "sec-7-9",
        "sec-7-1-3": "sec-7-10",
    },
    "08_speech-separation.md": {
        **{f"sec-8-1-{index}": f"sec-8-{index + 2}"
           for index in range(1, 7)},
    },
    "09_source-tracking.md": {
        "sec-u-f7a201edb9": "sec-9-6",
    },
    "10_engineering-practice.md": {
        "sec-u-38edc89a29": "sec-u-f9e36547b8",  # 有状态重采样旧标题
        "sec-u-a4719e76d5": "sec-u-7a1c8bc661",  # E10-02：两帧挂起 -> 两帧尾音保持
        "sec-10-9-1": "sec-10-8-1",              # 开源实现怎样用作工程对照
        "sec-10-12": "sec-10-11",                # 本章小结
    },
    "11_selection-guide.md": {
        "sec-u-65c140a10b": "sec-u-1cf1ad108e",  # E11-01旧题4标题
        "sec-u-1cf1ad108e": "sec-u-a08378b898",  # E11-01
        "sec-u-b2f4d38876": "sec-u-9348337294",  # E11-02
        "sec-u-994102bab0": "sec-u-c1060fa0f9",  # E11-03
        "sec-u-67b9b2107f": "sec-u-2499dabbff",  # E11-04
        "sec-u-b322947d6d": "sec-u-4284a6861e",  # E11-05
        "sec-u-7b9a84ce60": "sec-u-901b881949",  # E11-06
        "sec-u-a3b87b0c99": "sec-u-0ac305697e",  # E11-07
    },
    "12_appendix-symbols-math.md": {
        "sec-16-3-5": "sec-16-3",
        "sec-16-3-6": "sec-16-3",
        "sec-16-3-7": "sec-16-3",
        "sec-16-3-8": "sec-16-3-3",
        "sec-16-3-9": "sec-16-3-3",
        "sec-16-3-10": "sec-16-3-3",
        "sec-16-3-11": "sec-16-3-2",
        "sec-u-6ce2c569c7": "sec-16-4",          # old E16-07 pruned
        "sec-u-10a0b6ac2d": "sec-u-a2ed8e6072",  # old E16-08 -> E16-07
        "sec-u-e9ca3332f8": "sec-16-4",          # old E16-09 pruned
        "sec-u-723599188d": "sec-16-4",          # old E16-10 pruned
        "sec-u-d2443c8742": "sec-16-4",          # old E16-11 pruned
        "sec-u-058a6268ba": "sec-u-0b918e5d39",  # old E16-12 -> E16-08
        "sec-u-aebef34090": "sec-u-8fee69c59c",  # old E16-13 -> E16-09
        "sec-u-7bf1c7588b": "sec-u-1f33bf7dc5",  # old E16-14 -> E16-10
        "sec-u-a87d89ce6d": "sec-u-b7c6ebe600",  # old E16-15 merged into E16-14
        "sec-u-ff8521801e": "sec-u-ee4c3e8f61",  # old E16-16 -> E16-11
        "sec-u-7cb5badab0": "sec-u-76145f9ed0",  # old E16-17 -> E16-12
        "sec-u-401a933cc4": "sec-u-261df26518",  # old E16-18 -> E16-13
        "sec-u-0cf65b9988": "sec-u-b7c6ebe600",  # old E16-19 -> E16-14
        "sec-u-670d7ca59e": "sec-u-6f6ba9c2c8",  # old E16-20 -> E16-15
        "sec-u-0e76cc66ab": "sec-16-4",          # old E16-21 pruned
    },
    "13_appendix-guide.md": {
        "sec-u-b74f82c81d": "sec-u-36ad23ae3f",
        "sec-13-3-1": "sec-u-36ad23ae3f",
        "sec-u-69adc71a3e": "sec-u-7db9df573f",  # old E17-05 merged into E17-12
        "sec-u-ddd1562ee5": "sec-u-459edf7e20",  # old E17-06 -> E17-05
        "sec-u-fd4c1380d6": "sec-u-032701b632",  # old E17-07 -> E17-06
        "sec-u-84c390b4b7": "sec-u-471e180f3e",  # old E17-08 -> E17-07
        "sec-u-290348b839": "sec-u-16608f4cde",  # old E17-09 -> E17-08
        "sec-u-5e2cfdc863": "sec-u-69a1412d8e",  # old E17-10 merged into E17-11
        "sec-u-852b50db5d": "sec-u-990d0e2aed",  # old E17-11 -> E17-09
        "sec-u-a41bfbbab5": "sec-u-fe30a89239",  # old E17-12 -> E17-10
        "sec-u-d44d703e08": "sec-u-7db9df573f",  # old E17-13 merged into E17-12
        "sec-u-a237f010d4": "sec-u-69a1412d8e",  # old E17-14 -> E17-11
        "sec-u-def2767fe6": "sec-u-7db9df573f",  # old E17-15 merged into E17-12
        "sec-u-3396687cce": "sec-u-e78e69cadc",  # old E17-16 -> E17-13
    },
    "04_source_reproduction.md": {
        "sec-u-3b6143dc0a": "sec-u-0f9fac89c1",
        "8-固定分布式源码静态合同原函数与完整算法": "sec-u-0f9fac89c1",
        "sec-u-ce14b272d7": "sec-u-282fb01ca8",
        "原matlab静态合同与五项独立控制": "sec-u-282fb01ca8",
    },
    "05_exercises_and_audio.md": {
        "sec-u-c3cb1966de": "sec-u-94bcdb96bb",
        "61-第16章空间声场与可穿戴阵列对照实验": "sec-u-94bcdb96bb",
        "sec-u-3d1374a482": "sec-u-426c78e4cd",
        "62-第17章条件目标提取与近场距离气泡对照实验": "sec-u-426c78e4cd",
        "461-e13-11平均不变与逐通道等变": "sec-46-1",
        "461-e15-11平均不变与逐通道等变": "sec-46-1",
        "462-e13-12共同幅度尺度与安全能量对数": "sec-46-2",
        "462-e15-12共同幅度尺度与安全能量对数": "sec-46-2",
        "463-e13-13先移时间原点再恢复斜率单位": "sec-46-3",
        "463-e15-13先移时间原点再恢复斜率单位": "sec-46-3",
        "464-e13-14相同rir能量比为何产生不同输出": "sec-46-4",
        "464-e15-14相同rir能量比为何产生不同输出": "sec-46-4",
        "466-e13-15跨过拟合阈值仍可能缺少衰减尾部": "sec-46-6",
        "466-e15-15跨过拟合阈值仍可能缺少衰减尾部": "sec-46-6",
        "471-十八题分别核对什么": "sec-47-1",
        "482-25题的输入答案与代码控制": "sec-48-2",
        "sec-u-461361f106": "sec-u-63fb28005f",  # §52 E03-18 -> E03-16
        "sec-u-ed41d4e99e": "sec-u-0697cc2c95",  # §53 E04-24 -> E04-20
        "sec-u-935768b998": "sec-u-09a225e053",  # §54 E04-25 -> E04-21
    },
}


def historical_aliases(source_name: str, current_id: str) -> tuple[str, ...]:
    original_name = historical_source(source_name)
    renames = dict(HISTORICAL_SECTION_IDS.get(original_name, {}))
    renames.update(numbering_heading_aliases(source_name))
    aliases = {}
    for mapping in (renames, HISTORICAL_SEQUENTIAL_IDS.get(original_name, {})):
        for old_id, target in mapping.items():
            # The immutable original sequential map may name a subsequently
            # renamed semantic anchor. Follow that rename without rewriting
            # the old map or assigning its number to a different topic.
            visited = set()
            while target in renames:
                if target in visited:
                    raise ValueError(f"历史标题别名循环：{source_name}:{target}")
                visited.add(target)
                target = renames[target]
            if target == current_id:
                aliases[old_id] = None
    return tuple(aliases)


def has_historical_sequential_aliases(source_name: str) -> bool:
    return historical_source(source_name) in HISTORICAL_SEQUENTIAL_IDS
