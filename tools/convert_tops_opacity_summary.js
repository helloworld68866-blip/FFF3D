#!/usr/bin/env node

const fs = require("fs");
const path = require("path");

function decodeHtml(text) {
  return text
    .replace(/&nbsp;/g, " ")
    .replace(/&amp;/g, "&")
    .replace(/&lt;/g, "<")
    .replace(/&gt;/g, ">")
    .replace(/&quot;/g, '"')
    .replace(/&#39;/g, "'");
}

function extractCodeText(html) {
  const match = html.match(/<code[^>]*>([\s\S]*?)<\/code>/i);
  if (!match) {
    throw new Error("TOPS summary HTML does not contain a <code> results block");
  }
  return decodeHtml(match[1])
    .replace(/<br\s*\/?\s*>/gi, "\n")
    .replace(/<[^>]+>/g, "")
    .replace(/\r/g, "");
}

function numbersFromLine(line) {
  const matches = line.match(/[+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?/g);
  return matches ? matches.map(Number) : [];
}

function finiteNumber(value, label) {
  if (!Number.isFinite(value)) {
    throw new Error(`Expected finite number for ${label}`);
  }
  return value;
}

function key2(a, b) {
  return `${Number(a).toExponential(12)}|${Number(b).toExponential(12)}`;
}

function csvEscape(value) {
  const s = String(value);
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

function writeCsv(filePath, header, rows) {
  const out = [header.join(",")];
  for (const row of rows) {
    out.push(row.map(csvEscape).join(","));
  }
  fs.writeFileSync(filePath, `${out.join("\n")}\n`, "utf8");
}

function parseGrid(lines, marker, expectedCount) {
  const markerIndex = lines.findIndex((line) => line.includes(marker));
  if (markerIndex < 0) {
    throw new Error(`Missing grid marker: ${marker}`);
  }
  const values = [];
  for (let i = markerIndex + 1; i < lines.length && values.length < expectedCount; ++i) {
    values.push(...numbersFromLine(lines[i]));
  }
  if (values.length !== expectedCount) {
    throw new Error(`Expected ${expectedCount} values for ${marker}, found ${values.length}`);
  }
  return values.map((v, i) => finiteNumber(v, `${marker}[${i}]`));
}

function parseWarnings(lines) {
  const warnings = [];
  const warningIndex = lines.findIndex((line) =>
    line.includes("Requested densities listed below went off table boundaries"),
  );
  if (warningIndex < 0) {
    return warnings;
  }
  const compositionIndex = lines.findIndex((line) => line.includes("Normalized composition"));
  for (let i = warningIndex + 1; i < compositionIndex; ++i) {
    const values = numbersFromLine(lines[i]);
    if (values.length === 3) {
      warnings.push({
        temperature_keV: values[0],
        density_requested_g_cm3: values[1],
        density_used_g_cm3: values[2],
      });
    }
  }
  return warnings;
}

function parseComposition(lines) {
  const compositionIndex = lines.findIndex((line) => line.includes("Normalized composition"));
  if (compositionIndex < 0) {
    return [];
  }
  const entries = [];
  for (let i = compositionIndex + 1; i < lines.length; ++i) {
    if (lines[i].includes("Temperature grid used")) {
      break;
    }
    const parts = lines[i].trim().split(/\s+/);
    if (
      parts.length === 5 &&
      Number.isFinite(Number(parts[0])) &&
      Number.isFinite(Number(parts[1])) &&
      Number.isFinite(Number(parts[2]))
    ) {
      entries.push({
        number_fraction: Number(parts[0]),
        mass_fraction: Number(parts[1]),
        atomic_number: Number(parts[2]),
        chemical_symbol: parts[3],
        material_id: parts[4],
      });
    }
  }
  return entries;
}

function parseMeanOpacities(lines, densityCount, warningByPoint) {
  const rows = [];
  const headerPattern =
    /Density\s+Ross opa\s+Planck opa\s+No\. Free\s+Av Sq Free\s+T=\s*([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)/;

  for (let i = 0; i < lines.length; ++i) {
    const header = lines[i].match(headerPattern);
    if (!header) {
      continue;
    }
    const temperature = Number(header[1]);
    let readRows = 0;
    for (let j = i + 1; j < lines.length && readRows < densityCount; ++j) {
      const values = numbersFromLine(lines[j]);
      if (values.length !== 5) {
        continue;
      }
      const densityRequested = values[0];
      const warning = warningByPoint.get(key2(temperature, densityRequested));
      rows.push({
        temperature_keV: temperature,
        density_requested_g_cm3: densityRequested,
        density_used_g_cm3: warning ? warning.density_used_g_cm3 : densityRequested,
        density_was_clipped: Boolean(warning),
        rosseland_cm2_g: values[1],
        planck_cm2_g: values[2],
        free_electrons: values[3],
        avg_sq_free_electrons: values[4],
      });
      ++readRows;
    }
    if (readRows !== densityCount) {
      throw new Error(`Mean opacity block at T=${temperature} has ${readRows} rows`);
    }
  }
  return rows;
}

function parseMultigroupOpacities(lines, photonCount, warningByPoint) {
  const rows = [];
  const headerPattern =
    /Energy\s+Ross mg\s+Planck mg\s+for T,\s*density\s*=\s*([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s+([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)/;

  for (let i = 0; i < lines.length; ++i) {
    const header = lines[i].match(headerPattern);
    if (!header) {
      continue;
    }
    const temperature = Number(header[1]);
    const densityRequested = Number(header[2]);
    const warning = warningByPoint.get(key2(temperature, densityRequested));
    let readRows = 0;
    for (let j = i + 1; j < lines.length && readRows < photonCount; ++j) {
      const values = numbersFromLine(lines[j]);
      if (values.length !== 3) {
        continue;
      }
      rows.push({
        temperature_keV: temperature,
        density_requested_g_cm3: densityRequested,
        density_used_g_cm3: warning ? warning.density_used_g_cm3 : densityRequested,
        density_was_clipped: Boolean(warning),
        photon_energy_keV: values[0],
        rosseland_multigroup_cm2_g: values[1],
        planck_multigroup_cm2_g: values[2],
      });
      ++readRows;
    }
    if (readRows !== photonCount) {
      throw new Error(
        `Multigroup opacity block at T=${temperature}, density=${densityRequested} has ${readRows} rows`,
      );
    }
  }
  return rows;
}

function main() {
  const inputPath = process.argv[2];
  const outputDir = process.argv[3];
  if (!inputPath || !outputDir) {
    console.error("Usage: node tools/convert_tops_opacity_summary.js <TOPS summary.html> <output-dir>");
    process.exit(2);
  }

  const html = fs.readFileSync(inputPath, "utf8");
  const text = extractCodeText(html);
  const lines = text.split("\n").map((line) => line.trimEnd());

  const countLine = lines.find((line) => line.includes("Number of T"));
  const countMatch = countLine?.match(/Number of T\s*=\s*(\d+)\s+Number of rho\s*=\s*(\d+)\s+Number of materials\s*=\s*(\d+)/);
  if (!countMatch) {
    throw new Error("Missing TOPS count line");
  }
  const temperatureCount = Number(countMatch[1]);
  const densityCount = Number(countMatch[2]);
  const materialCount = Number(countMatch[3]);

  const materialLine = lines.find((line) => line.includes("TOPS results for"));
  const materialMatch = materialLine?.match(/TOPS results for\s+(.+?)\s+on\s+(.+)$/);
  const material = materialMatch ? materialMatch[1].trim() : "unknown";
  const topsDate = materialMatch ? materialMatch[2].trim() : "unknown";

  const temperatures = parseGrid(lines, "Temperature grid used the following", temperatureCount);
  const densities = parseGrid(lines, "Density grid used the following", densityCount);
  const photonEnergies = parseGrid(lines, "Photon grid used the following", 32);
  const warnings = parseWarnings(lines);
  const warningByPoint = new Map(warnings.map((w) => [key2(w.temperature_keV, w.density_requested_g_cm3), w]));
  const composition = parseComposition(lines);
  const meanRows = parseMeanOpacities(lines, densityCount, warningByPoint);
  const multigroupRows = parseMultigroupOpacities(lines, photonEnergies.length, warningByPoint);

  const expectedMeanRows = temperatureCount * densityCount;
  const expectedMultigroupRows = expectedMeanRows * photonEnergies.length;
  if (meanRows.length !== expectedMeanRows) {
    throw new Error(`Expected ${expectedMeanRows} mean rows, found ${meanRows.length}`);
  }
  if (multigroupRows.length !== expectedMultigroupRows) {
    throw new Error(`Expected ${expectedMultigroupRows} multigroup rows, found ${multigroupRows.length}`);
  }

  fs.mkdirSync(outputDir, { recursive: true });

  writeCsv(
    path.join(outputDir, "density_clipping_warnings.csv"),
    ["temperature_keV", "density_requested_g_cm3", "density_used_g_cm3"],
    warnings.map((w) => [w.temperature_keV, w.density_requested_g_cm3, w.density_used_g_cm3]),
  );

  writeCsv(
    path.join(outputDir, "mean_opacities.csv"),
    [
      "temperature_keV",
      "density_requested_g_cm3",
      "density_used_g_cm3",
      "density_was_clipped",
      "rosseland_cm2_g",
      "planck_cm2_g",
      "free_electrons",
      "avg_sq_free_electrons",
    ],
    meanRows.map((r) => [
      r.temperature_keV,
      r.density_requested_g_cm3,
      r.density_used_g_cm3,
      r.density_was_clipped,
      r.rosseland_cm2_g,
      r.planck_cm2_g,
      r.free_electrons,
      r.avg_sq_free_electrons,
    ]),
  );

  writeCsv(
    path.join(outputDir, "multigroup_opacities.csv"),
    [
      "temperature_keV",
      "density_requested_g_cm3",
      "density_used_g_cm3",
      "density_was_clipped",
      "photon_energy_keV",
      "rosseland_multigroup_cm2_g",
      "planck_multigroup_cm2_g",
    ],
    multigroupRows.map((r) => [
      r.temperature_keV,
      r.density_requested_g_cm3,
      r.density_used_g_cm3,
      r.density_was_clipped,
      r.photon_energy_keV,
      r.rosseland_multigroup_cm2_g,
      r.planck_multigroup_cm2_g,
    ]),
  );

  const metadata = {
    source_html: path.resolve(inputPath),
    material,
    tops_result_date: topsDate,
    units: {
      opacity: "cm2_per_g",
      temperature: "keV",
      density: "g_per_cm3",
      photon_energy: "keV",
    },
    counts: {
      temperatures: temperatureCount,
      densities: densityCount,
      materials: materialCount,
      photon_grid_points: photonEnergies.length,
      mean_rows: meanRows.length,
      multigroup_rows: multigroupRows.length,
      density_clipping_warning_rows: warnings.length,
    },
    composition,
    temperature_grid_keV: temperatures,
    density_grid_requested_g_cm3: densities,
    photon_grid_keV: photonEnergies,
    clipping_warning: {
      present: warnings.length > 0,
      interpretation:
        "TOPS reports that these requested density points were outside table boundaries and were evaluated at density_used_g_cm3.",
    },
    tables: {
      mean_opacities_csv: "mean_opacities.csv",
      multigroup_opacities_csv: "multigroup_opacities.csv",
      density_clipping_warnings_csv: "density_clipping_warnings.csv",
    },
    notes: [
      "The TOPS page labels the material as DT. Its normalized composition row reports Z=1, Chem. Sym.=H, Mat ID=4525; keep the TOPS material label and material id together when consuming this table.",
      "Rows with density_was_clipped=true should not be used as clean interpolation support points unless the provider intentionally accepts TOPS boundary clipping.",
      "The photon grid is exported as photon_energy_keV because the TOPS summary labels the column Energy; do not reinterpret it as group edges without a separate TOPS contract.",
    ],
  };
  fs.writeFileSync(path.join(outputDir, "metadata.json"), `${JSON.stringify(metadata, null, 2)}\n`, "utf8");

  const readme = `# TOPS DT Opacity Table\n\n` +
    `Source: ${path.resolve(inputPath)}\n\n` +
    `Material label: ${material}\n\n` +
    `TOPS date: ${topsDate}\n\n` +
    `Units: opacities in cm^2/g, temperature in keV, density in g/cm^3, photon energy in keV.\n\n` +
    `Generated files:\n\n` +
    `- \`metadata.json\`: grids, units, composition, counts, and caveats.\n` +
    `- \`mean_opacities.csv\`: Rosseland/Planck mean opacity and free-electron columns for each T/rho point.\n` +
    `- \`multigroup_opacities.csv\`: photon-energy-indexed Rosseland/Planck multigroup opacity rows for each T/rho point.\n` +
    `- \`density_clipping_warnings.csv\`: TOPS off-table density warnings.\n\n` +
    `Important caveats:\n\n` +
    `- ${warnings.length} requested density points were clipped by TOPS; use \`density_was_clipped\` and \`density_used_g_cm3\` when consuming tables.\n` +
    `- The source labels the material as DT, but the normalized composition row is Z=1 / H / Mat ID=4525. Preserve the material label and material id together.\n` +
    `- The photon grid is exported as \`photon_energy_keV\`; this converter does not claim those values are P3 group edges.\n`;
  fs.writeFileSync(path.join(outputDir, "README.md"), readme, "utf8");

  console.log(JSON.stringify(metadata.counts, null, 2));
}

main();
