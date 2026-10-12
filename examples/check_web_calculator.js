#!/usr/bin/env node
/* Cross-check the web calculator (web/app.js) against the R reference values in
 * examples/worked_scoring_example.csv.
 *
 * For each worked case this feeds the RAW inputs (age in months, sex, BMI,
 * waist, height, ALT) through the calculator's own pipeline:
 *   cdcExtBmiz(ageMonths, sex, bmi) -> ext_bmiz   (CDC extended BMI-for-age z)
 *   WHtR = waist / height                          (ratio, not percent)
 *   ln_ALT = log(ALT)                              (natural log)
 *   scoreYHSS / scoreNoALT -> lp, p
 * and compares lp/p against the R reference (analysis/score_frozen_model.R,
 * full-precision frozen coefficients) with a tolerance of 1e-9.
 *
 * Run from the repository root:
 *   node examples/check_web_calculator.js
 */
"use strict";
const fs = require("fs");
const path = require("path");

const YHSS = require(path.join(__dirname, "..", "web", "app.js"));

function parseCSV(text) {
  const lines = text.trim().split(/\r?\n/);
  const header = lines[0].split(",").map(s => s.replace(/^"|"$/g, ""));
  return lines.slice(1).map(line => {
    const cells = line.split(",").map(s => s.replace(/^"|"$/g, ""));
    const row = {};
    header.forEach((h, i) => { row[h] = cells[i]; });
    return row;
  });
}

const csv = fs.readFileSync(path.join(__dirname, "worked_scoring_example.csv"), "utf8");
const cases = parseCSV(csv);

const TOL = 1e-9;
let maxDpY = 0, maxDpN = 0, maxDlpY = 0, maxDlpN = 0, maxDz = 0, fail = 0;

console.log("id  age_m sex  max|dp_yhss|  max|dp_noalt|  verdict");
for (const c of cases) {
  const ageMonths = parseInt(c.age_months, 10);      // months, integer
  const sexInt = parseInt(c.sex_int, 10);            // 1=male, 2=female
  const sexF = sexInt === 2 ? 1 : 0;                 // female indicator
  const bmi = parseFloat(c.bmi);
  const whtr = parseFloat(c.waist_cm) / parseFloat(c.height_cm);
  const lnAlt = Math.log(parseFloat(c.alt_u_l));

  const bz = YHSS.cdcExtBmiz(ageMonths, sexInt, bmi);
  const inp = { ageMonths: ageMonths, sexF: sexF, extBmiz: bz.z, whtr: whtr };
  const sY = YHSS.scoreYHSS(Object.assign({ lnAlt: lnAlt }, inp));
  const sN = YHSS.scoreNoALT(inp);

  const dZ = Math.abs(bz.z - parseFloat(c.ext_bmiz));
  const dlpY = Math.abs(sY.lp - parseFloat(c.lp_yhss));
  const dpY = Math.abs(sY.p - parseFloat(c.p_yhss));
  const dlpN = Math.abs(sN.lp - parseFloat(c.lp_noalt));
  const dpN = Math.abs(sN.p - parseFloat(c.p_noalt));
  maxDz = Math.max(maxDz, dZ); maxDlpY = Math.max(maxDlpY, dlpY);
  maxDpY = Math.max(maxDpY, dpY); maxDlpN = Math.max(maxDlpN, dlpN);
  maxDpN = Math.max(maxDpN, dpN);
  const ok = dZ < 1e-9 && dlpY < TOL && dpY < TOL && dlpN < TOL && dpN < TOL;
  if (!ok) fail++;
  console.log(`${c.id}   ${c.age_months}    ${sexInt}    ${dpY.toExponential(2)}      ${dpN.toExponential(2)}       ${ok ? "PASS" : "FAIL"}`);

  // classification at the locked thresholds must also match
  const clsY = sY.p >= YHSS.MODEL.t90 ? 1 : 0;
  const clsN = sN.p >= YHSS.MODEL_NOALT.t90 ? 1 : 0;
  if (clsY !== parseInt(c.yhss_screen_positive, 10) ||
      clsN !== parseInt(c.noalt_screen_positive, 10)) {
    console.log(`  classification MISMATCH in case ${c.id}`); fail++;
  }
}

console.log("\nmax |dz(ext_bmiz)| =", maxDz.toExponential(2));
console.log("max |dlp| Y-HSS =", maxDlpY.toExponential(2), " no-ALT =", maxDlpN.toExponential(2));
console.log("max |dp|  Y-HSS =", maxDpY.toExponential(2), " no-ALT =", maxDpN.toExponential(2));
console.log("t90 (Y-HSS) =", YHSS.MODEL.t90, " t90 (no-ALT) =", YHSS.MODEL_NOALT.t90);
if (fail === 0) { console.log("\nALL CHECKS PASS: web calculator == R reference."); }
else { console.log(`\n${fail} CHECK(S) FAILED`); process.exit(1); }
