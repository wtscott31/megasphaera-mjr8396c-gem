#!/usr/bin/env python3
"""Lightweight integrity check for the manuscript repository."""
from pathlib import Path
import csv
import json
import math
import sys

ROOT = Path(__file__).resolve().parents[1]
required = [
    "models/Megasphaera_sp_MJR8396C_filled_anaerobic.xml",
    "models/Bifidobacterium_adolescentis_filled_anaerobic.xml",
    "media/MJR8396C_gut_like_medium.tsv",
    "media/MJR8396C_carbon_rich_medium.tsv",
    "media/MJR8396C_AA_rich_medium.tsv",
    "tables/Table_3_growth_and_exchange_three_media.csv",
    "tables/Table_5_micom_exchange_fraction_0.95.csv",
    "quality_control/Megasphaera_sp_MJR8396C_filled_anaerobic_memote.html",
]
for n in range(1, 8):
    required.append(f"figures/manuscript/Figure_{n}.pdf")

missing = [p for p in required if not (ROOT / p).exists()]
if missing:
    print("Missing required repository files:")
    for p in missing:
        print(" -", p)
    sys.exit(1)

summary = json.loads((ROOT / "analysis_summary.json").read_text())
expected_summary = {
    "total_reactions": 1849,
    "gpr_associated_reactions": 1335,
    "non_gpr_total": 514,
    "blocked_gut_like_medium": 766,
}
for key, expected in expected_summary.items():
    if summary.get(key) != expected:
        raise SystemExit(f"Unexpected {key}: {summary.get(key)} != {expected}")

with (ROOT / "tables/Table_3_growth_and_exchange_three_media.csv").open(newline="") as fh:
    rows = list(csv.DictReader(fh))
by_name = {r.get("Condition") or r.get("condition") or r.get("medium"): r for r in rows}

def numeric(row, *names):
    for name in names:
        if name in row and row[name] != "":
            return float(row[name])
    raise KeyError(names)

checks = {
    "Gut-like": 0.7005,
    "Carbohydrate-rich": 1.2006,
    "Carbon-rich": 1.2006,
    "Amino-acid-rich": 0.8327,
    "AA-rich": 0.8327,
}
seen = 0
for label, expected in checks.items():
    if label in by_name:
        value = numeric(by_name[label], "Growth rate (h^-1)", "Growth", "growth", "growth_rate")
        if not math.isclose(value, expected, rel_tol=0, abs_tol=5e-4):
            raise SystemExit(f"Growth check failed for {label}: {value}")
        seen += 1
if seen < 3:
    raise SystemExit("Could not identify all three manuscript media in Table 3")

print("Repository integrity checks passed.")
