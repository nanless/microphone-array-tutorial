import {readFileSync, writeFileSync} from "node:fs";
import {runInNewContext} from "node:vm";

const file = name => new URL(name, import.meta.url);
const html = readFileSync(file("array-lab.html"), "utf8");
const model = html.match(/<script id="lab-model">([\s\S]*?)<\/script>/)?.[1];
if (!model) throw new Error("Embedded lab model not found");
const Lab = runInNewContext(model + "\nLab;");
const check = process.argv.includes("--check");
for (const mode of ["single", "left", "right"]) {
  const {output, fractions} = Lab.render(Lab.synthesize(false), mode);
  const wav = Buffer.from(Lab.wav(output));
  const dest = file(mode + ".wav");
  if (check) {
    if (!readFileSync(dest).equals(wav)) throw new Error("Published WAV differs: " + mode);
  } else {
    writeFileSync(dest, wav);
  }
  console.log(mode, "power retained (%):", fractions.map(fraction => fraction * 100));
}
console.log("left, low B, power retained (%):", Lab.render(Lab.synthesize(true), "left").fractions.map(x => x * 100));
