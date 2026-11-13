/** Bundle pinned Nightingale web components and their notices for offline use. */
import { createHash } from "node:crypto";
import { readFile, readdir, writeFile } from "node:fs/promises";
import { dirname, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { build } from "esbuild";

const here = dirname(fileURLToPath(import.meta.url));
const assets = resolve(here, "../src/pocket_workbench/assets");
const bundle = "nightingale-components.js";
const notices = "nightingale-NOTICES.txt";
const built = await build({
  absWorkingDir: here,
  entryPoints: ["register.js"],
  bundle: true,
  minify: true,
  format: "iife",
  target: "es2020",
  legalComments: "eof",
  outfile: join(assets, bundle),
  metafile: true,
});
const dependencies = new Set();
for (const input of Object.keys(built.metafile.inputs)) {
  const parts = input.split("node_modules/");
  if (parts.length < 2) continue;
  const tail = parts.at(-1).split("/");
  const name = tail[0].startsWith("@") ? tail.slice(0, 2).join("/") : tail[0];
  dependencies.add(parts.slice(0, -1).join("node_modules/") + "node_modules/" + name);
}
let licenceText = "Notices for the packages bundled in nightingale-components.js\n";
for (const dependency of [...dependencies].sort()) {
  const dir = resolve(here, dependency);
  const { name, version, license } = JSON.parse(await readFile(join(dir, "package.json"), "utf8"));
  licenceText += `\n${"-".repeat(72)}\n${name} ${version} · ${license || "see notice"}\n`;
  const names = (await readdir(dir)).filter(file => /^(licen[cs]e|copying|notice)(\.|$)/i.test(file));
  for (const name of names.sort()) licenceText += `\n${await readFile(join(dir, name), "utf8")}\n`;
  if (!names.length) throw new Error(`Missing packaged licence for ${name}`);
}
await writeFile(join(assets, notices), licenceText, "utf8");
const checksums = join(assets, "checksums.json");
const hashes = JSON.parse(await readFile(checksums, "utf8"));
for (const file of [bundle, notices]) {
  hashes[file] = createHash("sha256").update(await readFile(join(assets, file))).digest("hex");
}
await writeFile(checksums, JSON.stringify(Object.fromEntries(Object.entries(hashes).sort()), null, 2) + "\n");
