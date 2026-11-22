"""Render actual chapter passages in memory; never rewrite publication outputs."""
from html.parser import HTMLParser
import math
from pathlib import Path
import re
import unittest

from scripts.build_site import render


ROOT = Path(__file__).resolve().parents[1]


class TopLevelLists(HTMLParser):
    def __init__(self):
        super().__init__()
        self.lists = []
        self.depth = 0
        self.in_item = False

    def handle_starttag(self, tag, attrs):
        if tag in ("ol", "ul"):
            if self.depth == 0:
                self.lists.append([])
            self.depth += 1
        elif tag == "li" and self.depth == 1:
            self.lists[-1].append("")
            self.in_item = True

    def handle_endtag(self, tag):
        if tag == "li" and self.depth == 1:
            self.in_item = False
        elif tag in ("ol", "ul"):
            self.depth -= 1

    def handle_data(self, data):
        if self.in_item:
            self.lists[-1][-1] += data


def chapter_fragment(filename, start, end):
    path = ROOT / "chapters" / filename
    source = path.read_text(encoding="utf-8")
    return source.split(start, 1)[1].split(end, 1)[0], path


def rendered_lists(source, path):
    html, _ = render(source, path)
    parser = TopLevelLists()
    parser.feed(html)
    return html, parser.lists


class SpatialParagraphTest(unittest.TestCase):
    def test_spatial_aliasing_has_three_separate_steps(self):
        source, path = chapter_fragment(
            "02_basics-signal-model.md", "**为什么 ", "等号边界可直接代入核对"
        )
        _, lists = rendered_lists("**为什么 " + source, path)
        self.assertEqual([len(items) for items in lists], [3])
        for item, phrase in zip(lists[0], ("相邻两麦", "相位以", "取值范围")):
            self.assertIn(phrase, item)
        self.assertIn(r"\psi = 2\pi d\sin\theta/\lambda", lists[0][0])
        self.assertIn(r"\sin\theta-\sin\theta_0=\frac{\lambda}{d}", lists[0][1])
        self.assertIn(r"\lambda/d>2", lists[0][2])

    def test_beamforming_derivation_and_postfilters_are_separate_list_items(self):
        cases = (
            ("#### 阵因子、旁瓣与白噪声增益上限", "上述 $3.01", ("主瓣与栅瓣", "第一旁瓣电平", "白噪声增益硬上限")),
            ("### 5.7 后置滤波", "#### 5.7.1", ("Zelinski", "McCowan", "对数谱幅度")),
        )
        for start, end, labels in cases:
            with self.subTest(section=start):
                source, path = chapter_fragment("05_beamforming.md", start, end)
                html, lists = rendered_lists(source, path)
                self.assertEqual([len(items) for items in lists], [len(labels)])
                for item, label in zip(lists[0], labels):
                    self.assertIn(label, item)
                self.assertNotIn("<pre>", html)

    def test_mvdr_robust_routes_are_four_separate_items(self):
        source, path = chapter_fragment(
            "05_beamforming.md", "除固定对角加载外，四条经典加固路线如下：", "排障对照表"
        )
        _, lists = rendered_lists(source, path)
        self.assertEqual([len(items) for items in lists], [4])
        for item, phrase in zip(lists[0], ("SOCP", "特征空间投影", "空间谱积分协方差重构", "前后向空间平滑")):
            self.assertIn(phrase, item)
        self.assertIn(r"\hat{\mathbf R}_{nn}", lists[0][0])

    def test_kalman_five_steps_keep_equations_and_explanations_in_their_items(self):
        source, path = chapter_fragment(
            "09_source-tracking.md", "**卡尔曼滤波（KF）的五步递推**", "方位角跨越表示区间"
        )
        html, lists = rendered_lists("**卡尔曼滤波（KF）的五步递推**" + source, path)
        self.assertEqual([len(items) for items in lists], [5])
        labels = ("预测状态", "预测协方差", "计算卡尔曼增益", "更新状态", "更新协方差")
        for number, (item, label) in enumerate(zip(lists[0], labels), start=3):
            self.assertIn(label, item)
            equation = rf"\tag{{9-{number}}}"
            self.assertEqual(item.count(equation), 1)
            self.assertEqual(html.count(equation), 1)
            self.assertEqual(item.count(r"\tag{"), 1)
        self.assertIn("未建模运动的过程噪声", lists[0][1])
        self.assertIn("较小的观测方差会增大角度增益", lists[0][2])
        self.assertIn("第一行的差称为新息", lists[0][3])
        self.assertNotIn("<pre>", html)

    def test_missing_intro_separator_reproduces_the_original_failure(self):
        source, path = chapter_fragment(
            "02_basics-signal-model.md", "**为什么 ", "等号边界可直接代入核对"
        )
        broken = ("**为什么 " + source).replace("：\n\n1.", "：\n1.", 1)
        _, lists = rendered_lists(broken, path)
        self.assertEqual(lists, [])

    def test_decimal_prose_is_not_a_list(self):
        _, lists = rendered_lists("频率为 3.4 kHz。距离为 0.05 m。", ROOT / "chapters/example.md")
        self.assertEqual(lists, [])

    def test_mask_beamforming_keeps_three_open_source_checks_separate(self):
        source, path = chapter_fragment(
            "05_beamforming.md", "上线前务必核对三处常见接口语义差异：", "### 5.10 波束形成算法总对比"
        )
        html, lists = rendered_lists(source, path)
        self.assertEqual([len(items) for items in lists], [3])
        self.assertIn("BAN 尺度与 RTF", lists[0][0])
        self.assertIn("SDW-MWF", lists[0][1])
        self.assertIn("声学成像功率图", lists[0][2])
        self.assertNotIn("<pre>", html)

    def test_reviewed_topic_changes_are_separate_rendered_paragraphs(self):
        checks = [
            ("02_basics-signal-model.md", "周期 Hann 窗", "会在两端补半窗"),
            ("02_basics-signal-model.md", "在掩码加权批量或流式统计中", "需要加载时采用"),
            ("04_doa-estimation.md", "相干多径降秩与前后向空间平滑", "root-MUSIC 多项式求根"),
            ("04_doa-estimation.md", "输入特征从原始复谱走向物理对齐", "输出表示从双分支解耦"),
            ("04_doa-estimation.md", "DCASE2022 官方基线", "不同届 DCASE Task 3 的输出格式不可混用"),
            ("04_doa-estimation.md", "不同届 DCASE Task 3 的输出格式不可混用", "至于单靠音频估计"),
            ("04_doa-estimation.md", "结果与不可观测边界", "本书原创纯 NumPy 基线位于"),
            ("05_beamforming.md", "最小值控制递归平均（MCRA", "经过帧对齐"),
            ("05_beamforming.md", "以参考麦", "虚抬降噪分贝"),
        ]
        for filename, first, second in checks:
            path = ROOT / "chapters" / filename
            html, _ = render(path.read_text(encoding="utf-8"), path)
            paragraphs = re.findall(r"<p>(.*?)</p>", html, flags=re.S)
            with self.subTest(filename=filename, first=first):
                containing_first = [p for p in paragraphs if first in p]
                containing_second = [p for p in paragraphs if second in p]
                self.assertTrue(containing_first)
                self.assertTrue(containing_second)
                self.assertFalse(any(second in p for p in containing_first))

    def test_noise_reduction_formula_has_its_own_block(self):
        source, path = chapter_fragment(
            "05_beamforming.md", "**工程选型三台阶**", "#### 5.10.1 延伸阅读"
        )
        html, _ = render("**工程选型三台阶**" + source, path)
        paragraphs = re.findall(r"<p>(.*?)</p>", html, flags=re.S)
        def paragraph_index(needle):
            matches = [index for index, paragraph in enumerate(paragraphs) if needle in paragraph]
            self.assertEqual(len(matches), 1, f"Expected one paragraph containing {needle}")
            return matches[0]

        intro = paragraph_index("工程选型三台阶")
        nr_intro = paragraph_index("以参考麦")
        nr_formula = paragraph_index(r"\tag{5-28}")
        nr_explanation = paragraph_index("虚抬降噪分贝")
        self.assertRegex(html, r"<h4[^>]*>同输入算例：零陷、WNG 与目标损伤</h4>")
        example = paragraph_index("取四麦 ULA")
        sinr_formula = paragraph_index(r"\tag{5-29}")
        summary = paragraph_index("四组数字浓缩了全章核心权衡")

        self.assertEqual(
            [intro, nr_intro, nr_formula, nr_explanation, example, sinr_formula, summary],
            sorted({intro, nr_intro, nr_formula, nr_explanation, example, sinr_formula, summary}),
            "Baseline, formulas, example and summary must stay in distinct ordered paragraphs",
        )
        self.assertTrue(paragraphs[nr_formula].strip().startswith("$$"))
        self.assertTrue(paragraphs[nr_formula].strip().endswith("$$"))
        self.assertTrue(paragraphs[sinr_formula].strip().startswith("$$"))
        self.assertTrue(paragraphs[sinr_formula].strip().endswith("$$"))
        self.assertLess(html.index(r"\tag{5-29}"), html.index("<table>"))
        self.assertLess(html.index("</table>"), html.index("四组数字浓缩了全章核心权衡"))

    def test_music_and_nearfield_formula_structure(self):
        source, path = chapter_fragment("04_doa-estimation.md", "思想三步", "**正交性成立的三步证明**")
        html, lists = rendered_lists("**MUSIC思想三步" + source, path)
        self.assertEqual([len(items) for items in lists], [3])
        self.assertIn(r"\tag{4-15}", lists[0][2])
        self.assertNotIn("<pre>", html)

        nf_source, nf_path = chapter_fragment(
            "04_doa-estimation.md", "**锚点算法：角度 × 距离球面扫描**", "**结果与不可观测边界**"
        )
        nf_html, _ = render("**锚点算法：角度 × 距离球面扫描**" + nf_source, nf_path)
        self.assertIn(r"\tag{4-29}", nf_html)
        self.assertIn("分母为 $3$", nf_html)
        self.assertNotIn("<pre>", nf_html)

    def test_mwf_cost_conditions_and_rank_one_relation_have_separate_blocks(self):
        start = "#### 5.5.1 多通道维纳滤波与秩一目标模型"
        source, path = chapter_fragment(
            "05_beamforming.md", start, "#### 5.5.2 Frost 波束形成器：时域约束更新"
        )
        html, lists = rendered_lists(start + source, path)
        self.assertEqual(lists, [])
        paragraphs = re.findall(r"<p>(.*?)</p>", html, flags=re.S)
        # Locate technical paragraphs by their complete distinctive content,
        # rather than requiring a historical list layout that no longer exists.
        checks = (r"s_{\mathrm{ref}}", r"\tag{5-13}", r"\mu\to\infty",
                  r"\tag{5-14}", "目标协方差满秩")
        indices = []
        for needle in checks:
            found = [i for i, paragraph in enumerate(paragraphs) if needle in paragraph]
            self.assertEqual(len(found), 1, needle)
            indices.append(found[0])
        self.assertEqual(indices, sorted(set(indices)))
        for tag in (r"\tag{5-13}", r"\tag{5-14}"):
            formula = next(p.strip() for p in paragraphs if tag in p)
            self.assertTrue(formula.startswith("$$") and formula.endswith("$$"))
        self.assertIn(r"\frac{\phi_s q}{\mu+\phi_s q}", html)
        self.assertIn(r"a_r=1", paragraphs[indices[2]])
        self.assertNotIn("<pre>", html)

    def test_accdoa_external_and_book_azimuths_are_complementary(self):
        source = (ROOT / "chapters/04_doa-estimation.md").read_text(encoding="utf-8")
        self.assertIn(r"\operatorname{atan2}(0.25,0.433)\approx 30^\circ", source)
        self.assertIn(r"顺时针计角为 $60^\circ$", source)
        # The published three-decimal vector is rounded, not an exact 30-degree ray.
        x, y = .433, .25
        external = math.degrees(math.atan2(y, x))
        book = math.degrees(math.atan2(x, y))
        self.assertAlmostEqual(external + book, 90)
        self.assertLess(abs(external - 30), .001)
        self.assertLess(abs(book - 60), .001)
        self.assertLess(abs(math.hypot(x, y) - .5), .0001)


if __name__ == "__main__":
    unittest.main()
