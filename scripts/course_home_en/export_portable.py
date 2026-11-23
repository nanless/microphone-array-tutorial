#!/usr/bin/env python3
"""Export the offline English course introduction HTML from course_home_en sources."""
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
            path = BLOB + "chapters_en/00_overview.md"
        elif path.startswith("../assets/readme/lab/array-lab.html"):
            suffix = path[len("../assets/readme/lab/array-lab.html"):]
            path = PAGES + "assets/readme/lab/array-lab.html" + suffix
        elif path.startswith("../"):
            path = BLOB + path[3:]
        elif path == "research/index.html":
            path = BLOB + "codes/chapters/ch00/research_en/README.md"
        elif path == "research/03_industrial_deployment.html" or re.fullmatch(
            r"\d{2}_[a-z-]+\.html#sec-\d+-\d+", path
        ):
            path = PAGES + "site_en/" + path
        elif re.fullmatch(r"\d{2}_[a-z-]+\.html", path):
            path = BLOB + "chapters_en/" + path.removesuffix(".html") + ".md"
        elif path.startswith(PAGES + "site/") or path.startswith(PAGES + "site_en/"):
            pass
        elif path != REPOSITORY and path not in CASE_SOURCES and not path.startswith("#"):
            raise ValueError(f"Unmapped portable link: {path}")
        return f'href="{escape(path, quote=True)}"'

    template = re.sub(r'href="([^"]+)"', link, template)
    template = re.sub(r' data-lang-peer="[^\"]+"', "", template)
    template = template.replace('src="../assets/readme/cases/', f'src="{PAGES}assets/readme/cases/')
    root_open = '<div class="course" id="course-home" tabindex="-1">'
    if not template.startswith(root_open) or not template.endswith("</div>\n"):
        raise ValueError("Course home EN template must be wrapped in div#course-home and end with '</div>\\n'")
    template = (
        '<main class="course" id="course-home" tabindex="-1">'
        + template[len(root_open):-len("</div>\n")]
        + "</main>\n"
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Microphone Array Signal Processing · From Fundamentals to Industrial Practice</title>
<style>*{{box-sizing:border-box}}body{{margin:0;font-family:system-ui,sans-serif}}{css}</style>
</head>
<body>{template}
<script>{script}</script>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    template, css, script = ((SOURCE / name).read_text(encoding="utf-8")
                             for name in ("index.html", "home.css", "home.js"))
    target = DEST / "course-intro-en.html"
    data = portable(template, css, script).encode("utf-8")
    if args.check:
        if not target.is_file() or target.read_bytes() != data:
            raise SystemExit(f"Offline English course asset differs: {target}")
    else:
        target.write_bytes(data)
    print("Offline English course checked" if args.check else "Offline English course exported")


if __name__ == "__main__":
    main()
