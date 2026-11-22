"""Render real chapter passages in memory to protect list membership."""

from html.parser import HTMLParser
from pathlib import Path
import re
import unittest

from scripts import build_site


ROOT = Path(__file__).resolve().parents[1]


class ListParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.lists = []
        self.open_lists = []
        self.open_items = []
        self.open_paragraphs = []

    def handle_starttag(self, tag, attrs):
        if tag in ("ol", "ul"):
            record = {"tag": tag, "items": [], "depth": len(self.open_lists)}
            self.lists.append(record)
            self.open_lists.append(record)
        elif tag == "li":
            item = {"text": "", "paragraphs": 0, "paragraph_texts": []}
            self.open_lists[-1]["items"].append(item)
            self.open_items.append(item)
        elif tag == "p" and self.open_items:
            item = self.open_items[-1]
            item["paragraphs"] += 1
            item["paragraph_texts"].append("")
            self.open_paragraphs.append(item)

    def handle_endtag(self, tag):
        if tag in ("ol", "ul"):
            self.open_lists.pop()
        elif tag == "li":
            self.open_items.pop()
        elif tag == "p" and self.open_paragraphs:
            self.open_paragraphs.pop()

    def handle_data(self, data):
        for item in self.open_items:
            item["text"] += data
        if self.open_paragraphs:
            self.open_paragraphs[-1]["paragraph_texts"][-1] += data


class ParagraphParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.paragraphs = []
        self.in_paragraph = False

    def handle_starttag(self, tag, attrs):
        if tag == "p":
            self.paragraphs.append("")
            self.in_paragraph = True

    def handle_endtag(self, tag):
        if tag == "p":
            self.in_paragraph = False

    def handle_data(self, data):
        if self.in_paragraph:
            self.paragraphs[-1] += data


class EnhancementParagraphTest(unittest.TestCase):
    def passage(self, filename, start, end):
        path = ROOT / "chapters" / filename
        source = path.read_text(encoding="utf-8")
        begin = source.index(start)
        return source[begin:source.index(end, begin)], path

    def render_items(self, source, path, list_tag="ol"):
        html, _ = build_site.render(source, path)
        parser = ListParser()
        parser.feed(html)
        roots = [record for record in parser.lists if record['depth'] == 0]
        self.assertEqual(len(roots), 1, "Expected one contiguous top-level list")
        self.assertEqual(roots[0]["tag"], list_tag)
        return roots[0]["items"]

    def assert_labels(self, items, labels):
        self.assertEqual(len(items), len(labels))
        for item, label in zip(items, labels):
            self.assertTrue(item["text"].strip().startswith(label), item["text"])

    def test_reviewed_explanations_keep_distinct_paragraphs(self):
        cases = [
            ("06_aec.md", "#### 6.2.4 IPNLMS：", "**两抽头手算。**", "比例归一化最小均方", 7),
            ("06_aec.md", "第一行描述路径的漂移或突变", "逐行细读式(6-21)", "下式是便于复算", 3),
            ("08_speech-separation.md", "**AuxIVA（", "设统计帧数为", "这里的包络相关", 4),
            ("09_source-tracking.md", "PHD 强度 $D_t(\\vec{x})$", "下面的递推针对", "PHD 避免显式枚举", 2),
            ("09_source-tracking.md", "其中 $\\vec{\\delta}^{(i)}", "#### 9.3.3 可执行基线", "这些方法分别处理关联", 2),
            ("09_source-tracking.md", "阅读和复用开源追踪库时", "MHT 已有可阅读的受限参考", "把雷达或平面视觉教程代码迁移到麦克风阵列方向追踪时", 2),
        ]
        for filename, start, end, second_start, expected_paragraphs in cases:
            with self.subTest(filename=filename, start=start):
                source, path = self.passage(filename, start, end)
                html, _ = build_site.render(source, path)
                parser = ParagraphParser()
                parser.feed(html)
                self.assertEqual(len(parser.paragraphs), expected_paragraphs)
                self.assertNotIn(second_start, parser.paragraphs[0])
                self.assertTrue(parser.paragraphs[1].startswith(second_start))

    def test_aec_diagnostics_preserve_separate_topics(self):
        source, path = self.passage("06_aec.md", "#### 6.8.3 故障诊断与常见工程误区", "#### 6.8.4 代码与开源实现指引")
        html, _ = build_site.render(source, path)
        parser = ParagraphParser()
        parser.feed(html)
        starts = ["现场遇到“回声消不干净或双讲吞字”", "采样率偏移（SRO）的累积账本。",
                  "结合 D1～D8，工程联调中最需警惕"]
        self.assertEqual(len(parser.paragraphs), len(starts))
        for paragraph, start in zip(parser.paragraphs, starts):
            self.assertTrue(paragraph.startswith(start), paragraph)
            for other in starts:
                if other != start:
                    self.assertNotIn(other, paragraph)

    def test_device_selection_and_wpe_directions_lists(self):
        cases = [
            ("06_aec.md", "下面逐个子图拆读图26", "#### 6.4.1 自动双讲检测",
             ["图26（a）", "图26（b）", "图26（c）"],
             None),
            ("07_wpe-dereverberation.md", "WPE 的改进主要围绕四个限制展开", "参数设置应先",
             ["分布式阵列", "麦克风很多时", "短块和强噪声下", "系统级方法"],
             "神经网络可以为这些解析结构"),
        ]
        for filename, start, end, labels, continuation in cases:
            with self.subTest(filename=filename):
                source, path = self.passage(filename, start, end)
                items = self.render_items(source, path, list_tag="ul")
                self.assert_labels(items, labels)
                if continuation is not None:
                    self.assertTrue(all(continuation not in item["text"] for item in items))
                    html, _ = build_site.render(source, path)
                    self.assertLess(html.index("</ul>"), html.index(continuation))

    def test_dtd_five_methods_and_their_continuations(self):
        source, path = self.passage("06_aec.md", "**常见 DTD 方法**", "**从二值硬冻结走向连续软步长**")
        items = self.render_items(source, path)
        self.assert_labels(items, ["Geigel 幅度比较", "向量型 Benesty NCC", "频域 NCC", "相干函数法", "双滤波器法"])
        self.assertEqual(items[0]["paragraphs"], 3)
        self.assertIn("经典示例阈值", items[0]["text"])
        self.assertEqual(items[1]["paragraphs"], 7)
        self.assertIn(r"\tag{6-27}", items[1]["text"])
        self.assertIn("WebRTC AEC3", items[4]["text"])

    def test_wpe_five_steps_preserve_update_order(self):
        source, path = self.passage("07_wpe-dereverberation.md", "每个频点的一轮工程迭代", "边界检查包括")
        items = self.render_items(source, path)
        self.assertEqual(len(items), 5)
        self.assertIn("初始化", items[0]["text"])
        self.assertIn("有效帧集合", items[1]["text"])
        self.assertIn("加载方式（如式(7-9)）求系数", items[2]["text"])
        self.assertIn("由新系数计算残差", items[3]["text"])
        self.assertIn("迭代上限", items[4]["text"])
        self.assertIn("不能自动宣称原似然每轮下降", source)

    def test_three_aec_beamforming_orders(self):
        source, path = self.passage("06_aec.md", "#### 6.6.1 麦克风阵列与 AEC 的三种级联顺序", "**波束切换残差手算。**")
        items = self.render_items(source, path)
        self.assert_labels(items, ["AEC→BF", "BF→AEC", "联合优化"])
        self.assertIn("同频带窄带近似", items[1]["text"])
        self.assertIn(r"H_{\mathrm{eff}}", items[1]["text"])
        self.assertNotIn("联合优化", items[1]["text"])

    def test_gss_three_parts_and_equation_stays_in_first_item(self):
        source, path = self.passage("08_speech-separation.md", "**GSS（导引源分离", "**紧凑实现顺序**")
        items = self.render_items(source, path)
        self.assert_labels(items, ["cACGMM 空间聚类", "导引（Guided）", "级联"])
        paragraphs = items[0]["paragraph_texts"]
        anchors = (
            "cACGMM 空间聚类",
            r"\tag{8-18}",
            "正定的形状矩阵",
            "混合模型的后验概率",
        )
        positions = [next(i for i, paragraph in enumerate(paragraphs) if anchor in paragraph)
                     for anchor in anchors]
        self.assertEqual(positions, sorted(set(positions)), "The four explanation tasks need distinct ordered paragraphs")
        self.assertTrue(paragraphs[positions[1]].strip().startswith("$$"), "The cACG density needs its own math block")
        self.assertTrue(paragraphs[positions[1]].strip().endswith("$$"))
        self.assertIn("其中", paragraphs[positions[2]])
        self.assertIn("软时频掩码", paragraphs[positions[3]])
        for item in items[1:]:
            self.assertNotIn(r"\tag{8-18}", item["text"])

    def test_eight_research_directions_preserve_arraydps(self):
        source, path = self.passage("17_appendix-guide.md", "### 17.3", "#### 专题：神经网络前端")
        headings = re.findall(r"(?m)^#### (研究方向[一二三四五六七八]：[^\n]+)$", source)
        sections = re.split(r"(?m)^#### 研究方向[一二三四五六七八]：[^\n]+\n", source)[1:]
        labels = ["预训练模型", "阵列无关", "目标说话人提取", "生成式语音增强",
                  "TF-GridNet", "几何无关前端", "CHiME-9", "扩散先验"]
        self.assertEqual(len(headings), 8)
        self.assertEqual(len(sections), 8)
        for heading, label in zip(headings, labels):
            self.assertIn(label, heading)
        html, _ = build_site.render(source, path)
        self.assertEqual(html.count("<h4"), 8, "Each independent direction needs a navigation heading")
        self.assertIn("FlowSep", sections[3])
        self.assertIn("DiCoW", sections[6])
        self.assertIn("ArrayDPS", sections[7])
        self.assertIn("而非由注册声纹指定身份", sections[2])
        self.assertIn("最小复现", sections[7])

    def test_twelve_debugging_cases_and_continuations(self):
        source, path = self.passage("17_appendix-guide.md", "### 17.4", "### 17.5")
        items = self.render_items(source, path)
        self.assert_labels(items, [
            "MUSIC", "GCC-PHAT", "波束形成后", "超指向", "AEC", "SRP-PHAT", "多麦录音", "仿真很好",
            "STFT 逆变换后", "MEMS 差分阵", "PDM / TDM", "实时 AEC",
        ])
        for cause in ('几何不可辨识', '空间混叠', '多径反射'):
            self.assertIn(cause, items[5]['text'])
        for index, expected in [(0, "再检查源数"), (2, "噪声协方差既可"), (4, "若线性路径已收敛"), (7, "位置容差")]:
            self.assertEqual(items[index]["paragraphs"], 2)
            self.assertIn(expected, items[index]["text"])

    def test_missing_indent_cannot_pass_gss_membership_check(self):
        source, path = self.passage("08_speech-separation.md", "**GSS（导引源分离", "**紧凑实现顺序**")
        broken = source.replace("\n    ", "\n   ")
        html, _ = build_site.render(broken, path)
        parser = ListParser()
        parser.feed(html)
        self.assertNotIn(r"\tag{8-18}", parser.lists[0]["items"][0]["text"])


if __name__ == "__main__":
    unittest.main()
