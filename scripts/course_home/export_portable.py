"""Export the offline course introduction and GitHub diagram from the course sources."""
import argparse
from base64 import b64encode
from html import escape
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
SOURCE = Path(__file__).resolve().parent
DEST = ROOT / "assets/readme"
REPOSITORY = "https://github.com/nanless/microphone-array-tutorial"
BLOB = REPOSITORY + "/blob/main/"
AUDIO = ("single.wav", "left.wav", "right.wav")


def portable(template, css, script):
    def link(match):
        path = match[1]
        if path == "#main-content":
            path = BLOB + "chapters/00_overview.md"
        elif path.startswith("../"):
            path = BLOB + path[3:]
        elif path == "research/index.html":
            path = BLOB + "codes/chapters/ch00/research/README.md"
        elif path == "real_audio/README.md":
            path = BLOB + "codes/chapters/ch02/real_audio/README.md"
        elif re.fullmatch(r"\d{2}_[a-z-]+\.html", path):
            path = BLOB + "chapters/" + path.removesuffix(".html") + ".md"
        elif path != REPOSITORY and not path.startswith("#"):
            raise ValueError(f"Unmapped portable link: {path}")
        return f'href="{escape(path, quote=True)}"'

    for name in AUDIO:
        original = f'../assets/readme/lab/{name}'
        encoded = b64encode((DEST / "lab" / name).read_bytes()).decode("ascii")
        for attr in ("src", "href"):
            token = f'{attr}="{original}"'
            if template.count(token) != 1:
                raise ValueError(f"Unexpected portable audio layout: {token}")
            template = template.replace(token, f'{attr}="data:audio/wav;base64,{encoded}"' +
                                        (f' download="{name}"' if attr == "href" else ""))
    template = re.sub(r'href="(?!data:)([^"]+)"', link, template)
    template = template.replace('<div class="course" id="course-home" tabindex="-1">',
                                '<main class="course" id="course-home" tabindex="-1">', 1)
    template = template.removesuffix("</div>\n") + "</main>\n"
    return f"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>把声音听出方向 · 麦克风阵列入门</title>
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
.course-from-right,.course-wave,.course-hits circle{display:none}
.course-emitters circle{fill:#cbd7e7}
.course-emitters circle:first-child{fill:#ad402c}
.course-emitters text{fill:#17273d}
.course-emitters text:nth-child(3){fill:white}
"""
    if not still:
        style += """
@media (prefers-reduced-motion:no-preference) {
.course-wave{display:block;opacity:0;transform-origin:72px 144px;animation:spread 2.6s linear 3}
.course-hits circle{display:block;opacity:0;animation:hit 2.6s ease-out 3}
.course-hits circle:nth-child(1){animation-delay:.45s}
.course-hits circle:nth-child(2){animation-delay:.84s}
.course-hits circle:nth-child(3){animation-delay:1.22s}
.course-hits circle:nth-child(4){animation-delay:1.61s}
@keyframes spread{0%{opacity:0;transform:scale(.22)}12%{opacity:.8}65%{opacity:.62}80%,100%{opacity:0;transform:scale(2.35)}}
@keyframes hit{0%,19%,100%{opacity:0}5%{opacity:1}}
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
