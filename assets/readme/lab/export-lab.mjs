import {readFileSync, writeFileSync} from "node:fs";
import {runInNewContext} from "node:vm";

const file = name => new URL(name, import.meta.url);
const html = readFileSync(file("array-lab.html"), "utf8");
const model = html.match(/<script id="lab-model">([\s\S]*?)<\/script>/)?.[1];
if (!model) throw new Error("Embedded lab model not found");
const Lab = runInNewContext(model + "\nLab;");
const check = process.argv.includes("--check");
function exportAsset(name, bytes) {
  const dest = file(name);
  if (check) {
    if (!readFileSync(dest).equals(bytes)) throw new Error("Published asset differs: " + name);
  } else {
    writeFileSync(dest, bytes);
  }
}
const hero = readFileSync(file("../listening-room.svg"), "utf8");
const still = hero.replace(/<style>[\s\S]*?<\/style>/, "<style>.pulse { opacity: 0; }</style>");
if (still === hero) throw new Error("Hero animation stylesheet not found");
exportAsset("../listening-room-still.svg", Buffer.from(still));
for (const mode of ["single", "left", "right"]) {
  const {output, fractions} = Lab.render(Lab.synthesize(false), mode);
  exportAsset(mode + ".wav", Buffer.from(Lab.wav(output)));
  console.log(mode, "power retained (%):", fractions.map(fraction => fraction * 100));
}
console.log("left, low B, power retained (%):", Lab.render(Lab.synthesize(true), "left").fractions.map(x => x * 100));
