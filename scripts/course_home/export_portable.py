"""Export the offline course introduction and GitHub diagram from the course sources."""
import argparse
from html import escape
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
SOURCE = Path(__file__).resolve().parent
DEST = ROOT / "assets/readme"
REPOSITORY = "https://github.com/nanless/microphone-array-tutorial"
BLOB = REPOSITORY + "/blob/main/"
PAGES = "https://nanless.github.io/microphone-array-tutorial/"
CASE_SOURCES = {
    "https://zenodo.org/records/3630471",
    "https://arxiv.org/pdf/1909.01008#page=14",
    "https://doi.org/10.1109/MSP.2012.2205285",
}


def portable(template, css, script):
    def link(match):
        path = match[1]
        if path == "#main-content":
            path = BLOB + "chapters/00_overview.md"
        elif path == "../assets/readme/lab/array-lab.html":
            path = PAGES + "assets/readme/lab/array-lab.html"
        elif path.startswith("../"):
            path = BLOB + path[3:]
        elif path == "research/index.html":
            path = BLOB + "codes/chapters/ch00/research/README.md"
        elif path == "real_audio/README.md":
            path = BLOB + "codes/chapters/ch02/real_audio/README.md"
        elif path == "research/03_industrial_deployment.html" or re.fullmatch(
            r"\d{2}_[a-z-]+\.html#sec-\d+-\d+", path
        ):
            path = PAGES + "site/" + path
        elif re.fullmatch(r"\d{2}_[a-z-]+\.html", path):
            path = BLOB + "chapters/" + path.removesuffix(".html") + ".md"
        elif path != REPOSITORY and path not in CASE_SOURCES and not path.startswith("#"):
            raise ValueError(f"Unmapped portable link: {path}")
        return f'href="{escape(path, quote=True)}"'

    template = re.sub(r'href="([^"]+)"', link, template)
    template = template.replace('src="../assets/readme/cases/', f'src="{PAGES}assets/readme/cases/')
    root_open = '<div class="course" id="course-home" tabindex="-1">'
    if not template.startswith(root_open) or not template.endswith("</div>\n"):
        raise ValueError("Course home template must be wrapped in div#course-home and end with '</div>\\n'")
    template = (
        '<main class="course" id="course-home" tabindex="-1">'
        + template[len(root_open):-len("</div>\n")]
        + "</main>\n"
    )
    return f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>麦克风阵列信号处理 · 从入门到工业实践</title>
<style>*{{box-sizing:border-box}}body{{margin:0;font-family:system-ui,sans-serif}}{css}</style>
</head>
<body>{template}
<script>{script}</script>
</body>
</html>
"""


def cover(template, still):
    svg = re.search(r"<svg\b[\s\S]*?</svg>", template).group()
    svg = svg.replace("<svg ", '<svg xmlns="http://www.w3.org/2000/svg" ', 1)
    style = """
.course-out-b,.course-out-caption-right,.course-polar-left,.course-polar-right,.course-polar-markers-left,.course-polar-markers-right,.course-wave,.course-hits circle{display:none}
"""
    if not still:
        style += """
@media (prefers-reduced-motion:no-preference) {
.course-wave-a{display:block;opacity:0;transform-origin:104px 108px;animation:spread 2.8s linear 3}
.course-hits circle{display:block;opacity:0;animation:hit 2.8s ease-out 3}
.course-hits circle:nth-child(1){animation-delay:.40s}
.course-hits circle:nth-child(2){animation-delay:.76s}
.course-hits circle:nth-child(3){animation-delay:1.12s}
.course-hits circle:nth-child(4){animation-delay:1.48s}
@keyframes spread{0%{opacity:0;transform:scale(.35)}14%{opacity:.92}70%{opacity:.58}84%,100%{opacity:0;transform:scale(2.45)}}
@keyframes hit{0%,20%,100%{opacity:0}6%{opacity:1}}
}
"""
    return svg.replace("</defs>", "</defs>\n<style>" + style + "</style>", 1) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    template, css, script = ((SOURCE / name).read_text(encoding="utf-8")
                             for name in ("index.html", "home.css", "home.js"))
    for name, output in {
        "course-intro.html": portable(template, css, script),
        "course-cover.svg": cover(template, still=False),
        "course-cover-still.svg": cover(template, still=True),
    }.items():
        target = DEST / name
        data = output.encode("utf-8")
        if args.check:
            if target.read_bytes() != data:
                raise SystemExit(f"Offline course asset differs: {target}")
        else:
            target.write_bytes(data)
    print("Offline course and GitHub diagram checked" if args.check else "Offline course and GitHub diagram exported")


if __name__ == "__main__":
    main()
