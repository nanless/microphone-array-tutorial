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
    },
    "15_distributed-enhancement.md": {
        "sec-u-d261d70382": "sec-u-89268d58ca",  # 作者MATLAB静态合同 -> 工程实现细节
        "sec-u-b345641c53": "sec-u-a6af0cfe91",  # 可听控制 -> 可听对照音频
        "sec-u-320ff7d3af": "sec-u-87bec6abd7",  # 采样率偏差控制 -> 采样率偏差补偿
        "sec-u-ae60fe614c": "sec-u-c8f4d1ccc4",  # E13-01 全角空格 -> 冒号
        "sec-u-d2277f09c0": "sec-u-b01a879a4a",  # E13-02
        "sec-u-990fd59df5": "sec-u-db8a2119c9",  # E13-03
        "sec-u-9493aad880": "sec-u-ada0f6befb",  # E13-04
        "sec-u-a123bcf090": "sec-u-275e8f69e6",  # E13-05
        "sec-u-b7dac779da": "sec-u-0f09bd925b",  # E13-06
        "sec-u-9bcc32b7fe": "sec-u-66bc517dc7",  # E13-07
        "sec-u-40ecac23db": "sec-u-86e56e7136",  # E13-08
        "sec-u-15a9194ba1": "sec-u-2a988d60f2",  # E13-09
        "sec-u-ecaa4b04e4": "sec-u-76c29b9e9a",  # E13-10
        "sec-u-5159378837": "sec-u-8f41df6077",  # E13-11
        "sec-u-481fa3de3e": "sec-u-1ba11b3526",  # E13-12
        "sec-u-43c36175b8": "sec-u-55891f9c74",  # E13-13
        "sec-u-478a33071c": "sec-u-f54078f6df",  # E13-14
        "sec-u-4bad0eb87f": "sec-u-7cead60afd",  # E13-15
        "sec-u-d6bfdb5b5e": "sec-u-9f6977a658",  # E13-16
        "sec-u-270f9db859": "sec-u-2ab25b15a4",  # E13-17
        "sec-u-4e97b65795": "sec-u-03e1ca65d2",  # E13-18
        "sec-u-0c9a031714": "sec-u-87f37a1d30",  # E13-19
        "sec-u-5bc0cc2c13": "sec-u-0dcd9a42fe",  # E13-20
        "sec-u-9e15aaffec": "sec-u-2b170d53fc",  # E13-21
        "sec-u-4db262b2cd": "sec-u-0469bde73b",  # E13-22
        "sec-u-0b8e41071b": "sec-u-c7f8fcc136",  # E13-23
        "sec-u-7e45f79bb0": "sec-u-f7034ff89c",  # E13-24
        "sec-u-632be6e3e5": "sec-u-13a2e4109a",  # E13-25
        "sec-u-60abc8d82d": "sec-u-1979c49802",  # E13-26
        "sec-u-1f2fb8b5ba": "sec-u-8f5d54d722",  # E13-27
    },
    "00_overview.md": {
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
        "sec-u-05ec932a34": "sec-u-c36b7ac16e",  # 客厅里的远场语音问题：标题改为远距离
        "sec-1-1": "sec-u-1d5f9c4fa9",           # 1.1 双耳听觉的三类定位线索 -> 统一无编号节
    },
    "04_doa-estimation.md": {
        "sec-u-442f837dd2": "sec-u-5d3aebe0fb",  # MUSIC 旧标题
        "sec-u-8dda6571c9": "sec-u-e9bdae6530",  # ESPRIT 旧标题
        "sec-u-1e0973b8ca": "sec-u-eb21d54fe3",  # 算例 4-3 R -> R_{xx}
        "sec-4-6-1": "sec-u-08146bdaaf",  # 原 4.6.1 宽带聚焦
        "sec-4-7-1": "sec-u-00305faff5",  # 近场模型失配
        "sec-4-7-2": "sec-u-374bae76e9",  # 可执行基线与外部实现
    },
    "05_beamforming.md": {
        "sec-u-fb93744989": "sec-5-5-1",  # 多通道维纳滤波与秩一目标模型 -> 5.5.1
        "sec-u-b11411e989": "sec-5-5-2",  # Frost 波束形成器：时域约束更新 -> 5.5.2
    },
    "03_array-geometry.md": {
        "sec-u-5834a7d67e": "sec-3-3-1",
        "sec-u-20a0d88f74": "sec-3-3-2",
        "sec-u-d14b5899a5": "sec-3-3-3",
    },
    "06_aec.md": {
        "sec-6-1-3": "sec-6-2",
        **{f"sec-6-1-{index}": f"sec-6-{index - 1}"
           for index in range(4, 19)},
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
    },
    "11_selection-guide.md": {
        "sec-u-65c140a10b": "sec-u-1cf1ad108e",  # E11-01旧题4标题
    },
    "13_appendix-guide.md": {
        "sec-u-b74f82c81d": "sec-u-36ad23ae3f",
        "sec-13-3-1": "sec-u-36ad23ae3f",
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
