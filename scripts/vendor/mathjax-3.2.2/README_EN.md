# MathJax 3.2.2 and SRE 4.0.6: Offline Formula Assets for the Combined PDF

**[简体中文](README.md)**

This directory supplies offline math-rendering assets for `scripts/build_pdf.py` when printing the combined monograph PDF. Formulas are rendered into SVG glyph paths by MathJax 3.2.2 and annotated with complete MathML structural descriptions using the English rule set of Speech Rule Engine (SRE) 4.0.6. The static course website continues to load its online formula assets through its own site builder. Versions here are pinned strictly for deterministic build reproduction; we do not claim that they are the latest upstream releases.

## Official Provenance and Licenses

MathJax is vendored from the [official npm tarball](https://registry.npmjs.org/mathjax/-/mathjax-3.2.2.tgz) with archive SHA-256 `1b9c0a1c44df864e915690558e72adb9cc5203360daefd385084ced3b6c64c09`; its verbatim Apache-2.0 license text is preserved in [LICENSE](LICENSE). All asset files retain their exact upstream bytes, stripping only the leading `es5/` directory prefix from the archive paths.

The English speech rules come from the [official SRE 4.0.6 npm tarball](https://registry.npmjs.org/speech-rule-engine/-/speech-rule-engine-4.0.6.tgz), with archive SHA-256 `a8ca04e68178cd934ac57866c30af78c246fcba448ef590071398ac6486e0846`. MathJax's `a11y/sre.js` embeds SRE 4.0.6, and the bundled rule maps match that exact version. The verbatim Apache-2.0 license text for SRE is stored separately in [SRE-LICENSE](SRE-LICENSE) and is byte-identical to the [license at the pinned official commit](https://github.com/Speech-Rule-Engine/speech-rule-engine/blob/ba19057d31620aa1f7985893ad5cff090bc70e6d/LICENSE). No trained models or external browser runtime dependencies have been introduced.

| Local File | Role | SHA-256 |
|---|---|---|
| `tex-svg.js` | TeX parser and SVG glyph path renderer | `d4295dc33744836935c1399feece5159577b34c5c8ffb9f1c6324cd82e03a882` |
| `input/tex/extensions/boldsymbol.js` | Bold math symbol extension used across the book | `d6771fee0772db2657796c8d0e20e1878bb3237f6d3ed1e828e1834a4ff743ca` |
| `a11y/sre.js` | Browser interface embedding SRE 4.0.6 | `a41dc487037b75d6dae75677e36fe2f03d1d2684f232d96f0422a11b35051660` |
| `sre/mathmaps/base.json` | Base structural rule map | `83311887b069476a6a4d0c6afe7f8e55075c85b44b029588739c0d835fc912b8` |
| `sre/mathmaps/en.json` | English speech rule map | `19d9b309dc0d25d2bde2d3130bd3772b6d7c3f8316ebbd85306bef55ac89b2a7` |
| `LICENSE` | MathJax Apache-2.0 license | `cfc7749b96f63bd31c3c42b5c471bf756814053e847c10f3eb003417bc523d30` |
| `SRE-LICENSE` | Speech Rule Engine Apache-2.0 license | `0d542e0c8804e39aa7f37eb00da5a762149dc682d7829451287e11b938e94594` |
| `output/chtml/fonts/woff-v2/MathJax_Math-Regular.woff` | SVG layout unit metrics (does not draw SVG glyph paths) | `c01d3321e89b403c4b811aa153c4e618eda3421f92d8a072a02c8d190782a191` |

## Build Verification and Assistive-Text Boundaries

Run `.venv/bin/python scripts/build_pdf.py --build-date 2026-10-05` from the repository root; see the [Build Toolchain Guide](../../README_EN.md) for the complete workflow. `check_mathjax_assets` verifies the full SHA-256 digests of the main SVG script, the `boldsymbol` extension, the SRE browser interface, rule maps, SRE license, and the unit-metric font, and additionally checks the presence of the MathJax license file. Any missing file or digest mismatch aborts the build immediately. When headless Chrome loads the metric font and the `base`/`en` rule maps at runtime, it re-verifies their full byte digests and tests rule initialization against a fixed fraction control—merely awaiting `ready` or checking for a non-empty string would miss silent degradation when `base.json` is absent.

SVG path dimensions are expressed in TeX `ex` units, whose nominal x-height is `0.442em`. Inheriting CJK body fonts directly would inflate formula metrics and trigger whole-book page scaling during printing; therefore, the original `MathJax_Math-Regular.woff` is applied strictly to formula pre-measurement and SVG containers, with synthetic bold/italic disabled and actual font loading awaited prior to layout. Body prose remains `16px`, formulas inherit the font size of their enclosing paragraph or table cell, and CJK `<text>` nodes inside SVGs use the CJK body font. This WOFF file serves solely for unit and width calibration—visible mathematical glyphs are still drawn from SVG vector paths, and forcibly shrinking outer SVG bounding boxes to fit page width is never treated as valid mathematical typesetting.

Every formula occurrence is bound to its exact TeX source digest. Inside the browser, complete MathML is exported from the live MathJax math tree to synthesize a non-empty English structural description, and accessibility attributes are attached exclusively to real connected SVG elements. After PDF printing, the build verifies every formula's source identity against a unique real MCID and non-empty visible drawing operations, preserves the parent tree and glyphs, rewrites the corresponding structure element from `Figure` to `Formula`, and writes `ActualText` onto both the structure element and the matching marked-content sequence. Missing formulas, duplicate identities, undefined TeX macros, empty descriptions, or unrendered formulas all fail release verification; formulas are never marked as decorative artifacts to bypass checks.

SRE 4.0.6 does not ship Chinese speech rules, so generated formula descriptions are tagged `en`. Symbols such as asterisks, overbars, and `H` superscripts are described according to their visible notation, whereas domain semantics (such as complex conjugation or Hermitian transpose) remain defined by the surrounding prose. Source formulas containing Chinese text labels in the Chinese edition may still appear verbatim inside English structural descriptions, so mixed-language pronunciation and alignment with prose notation definitions require separate human verification. Machine-generated structural descriptions do not prove complete reading order, domain semantics, or end-user assistive-technology usability; automated PDF/UA validation and manual screen-reader testing must be recorded separately.

The PDF builder marks untagged typographic space glyphs as `Layout` artifacts only when jointly proven by an embedded TrueType font, an `Identity` CID-to-GID mapping, a unique `ToUnicode` mapping to `U+0020`, zero glyph contours, and a non-zero horizontal advance—preserving original text positioning without ever reclassifying visible formula paths or unknown content as decorative artifacts. Details of this bounded post-processing and page-content presence checks appear in the [Build Toolchain Guide](../../README_EN.md) and do not constitute a blanket PDF/UA conformance claim.

## Retained Historical Assets

The legacy `tex-mml-chtml.js` bundle and 23 `output/chtml/fonts/woff-v2/MathJax_*.woff` files are retained byte-for-byte for historical build reproduction. The legacy main script SHA-256 is `300480069078b5892d2363a2b65e2dfbbf30fe5c80f83edbfecf4610fd093862`. One of those original WOFF files (`MathJax_Math-Regular.woff`) now also provides the SVG unit metrics described above, while the remaining 22 WOFF files are kept for historical CHTML reproduction. Compatibility checks verify the presence and non-empty size of all 23 WOFF files, while verifying the full SHA-256 digest specifically for the active metric font; we do not claim byte-level digest verification for the other 22 legacy WOFF files.

The PDF build's input digest hashes every file in this directory alongside the formula helper scripts so any change to vendored inputs invalidates stale builds; however, a source digest alone does not prove upstream provenance or license compliance. Whenever versions or file selections change, independently verify the official upstream archives, licenses, and rendered output, and never hand-edit minified scripts, rule maps, or font binaries.
