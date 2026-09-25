#!/usr/bin/env node
/* Check that the scenarios hardcoded in web/ventsim.html mean the same thing as
 * the Python ones.  The HTML leaves out parameters that are equal to a default,
 * so the comparison is done on *resolved* patients: build the patient with the
 * JS model and compare every field against the Python reference JSON.
 *
 *   python3 tools/dump_scenarios.py /tmp/scen.json && node tools/check_scenarios.js /tmp/scen.json
 */
const fs = require("fs");
const path = require("path");

const V = require(path.join(__dirname, "..", "web", "ventsim.js"));

const htmlPath = path.join(__dirname, "..", "web", "ventsim.html");
const refPath = process.argv[2] || "/tmp/scen.json";

function extractArray(text, marker) {
  const start = text.indexOf(marker);
  if (start < 0) throw new Error(`marker not found: ${marker}`);
  const i = text.indexOf("[", start);
  let depth = 0, quote = null, esc = false;
  for (let j = i; j < text.length; j++) {
    const c = text[j];
    if (quote) {
      if (esc) { esc = false; continue; }
      if (c === "\\") { esc = true; continue; }
      if (c === quote) quote = null;
      continue;
    }
    if (c === "'" || c === '"' || c === "`") { quote = c; continue; }
    if (c === "[") depth++;
    else if (c === "]") { depth--; if (depth === 0) return text.slice(i, j + 1); }
  }
  throw new Error("unbalanced array");
}

const html = fs.readFileSync(htmlPath, "utf8");
const scenarios = eval(extractArray(html, "var SCENARIOS = ["));
const ref = JSON.parse(fs.readFileSync(refPath, "utf8"));

const js = {};
for (const s of scenarios) js[s.key] = s;

let bad = 0;
const want = Object.keys(ref).sort();
const have = Object.keys(js).sort();
if (want.join(",") !== have.join(",")) {
  console.log("KEY MISMATCH\n  python:", want.join(","), "\n  js    :", have.join(","));
  bad++;
}

for (const key of want) {
  const a = ref[key], b = js[key];
  if (!b) { console.log(`${key}: missing from HTML`); bad++; continue; }
  if (a.title !== b.title) console.log(`${key}: title "${a.title}" != "${b.title}"`), bad++;
  if (a.adapt_hours !== (b.adapt_hours || 0))
    console.log(`${key}: adapt_hours ${a.adapt_hours} != ${b.adapt_hours}`), bad++;
  if (b.story !== undefined && a.story !== b.story) console.log(`${key}: story differs`), bad++;
  if (b.expect !== undefined && a.expect !== b.expect) console.log(`${key}: expect differs`), bad++;

  const pjs = V.makePatient(b.patient);
  for (const f of Object.keys(a.patient)) {
    const av = a.patient[f], jv = pjs[f];
    if (f === "lung") continue;
    if (av === null) {
      if (jv !== null && jv !== undefined) console.log(`${key}: patient.${f} null != ${jv}`), bad++;
      continue;
    }
    if (typeof av === "number") {
      if (typeof jv !== "number" || Math.abs(av - jv) > 1e-9)
        console.log(`${key}: patient.${f} ${av} != ${jv}`), bad++;
    } else if (typeof av === "object") {
      for (const kk of Object.keys(av)) {
        const jv2 = jv ? jv[kk] : undefined;
        if (typeof jv2 !== "number" || Math.abs(av[kk] - jv2) > 1e-9)
          console.log(`${key}: patient.${f}.${kk} ${av[kk]} != ${jv2}`), bad++;
      }
    } else if (f !== "notes" && String(av) !== String(jv)) {
      console.log(`${key}: patient.${f} ${JSON.stringify(av)} != ${JSON.stringify(jv)}`), bad++;
    }
  }

  const vjs = V.makeVent(b.vent);
  for (const vk of ["mode", "fio2", "peep", "vt_ml_kg", "rr", "ie", "dp"]) {
    const av = a.vent[vk], jv = vjs[vk];
    if (vk === "mode") { if (String(av) !== String(jv)) console.log(`${key}: vent.mode ${av} != ${jv}`), bad++; }
    else if (typeof jv !== "number" || Math.abs(av - jv) > 1e-9)
      console.log(`${key}: vent.${vk} ${av} != ${jv}`), bad++;
  }

  // and the operating point itself, at the scenario's own settings
  const res = V.solveSteady(pjs, vjs, { hco3: a.hco3 });
  void res;
}
console.log(bad === 0 ? `scenarios agree (${want.length} cases)` : `${bad} mismatch(es)`);
process.exit(bad ? 1 : 0);
