#!/usr/bin/env python3
"""Dependency-free static verification of the published v6 repository structure and supplemental file map."""
from pathlib import Path
import hashlib,csv,json,sys
R=Path(__file__).resolve().parents[1]
index=R/'supplementary/SUPPLEMENTARY_INDEX.csv'
assert index.exists(),index
rows=list(csv.DictReader(index.open()))
assert len(rows)==87,f'Expected 87 numbered source tables, got {len(rows)}'
missing=[]
for row in rows:
 for fn in row['canonical_files'].split(';'):
  if fn and not (R/fn).is_file():missing.append(fn)
for n in range(1,8):
 f=R/f'figures/manuscript/Figure_{n}.pdf'
 if not f.is_file():missing.append(str(f))
for f in ['models/original/Megasphaera_sp_MJR8396C_filled_anaerobic.xml',
          'models/MJR8396C_matched_sequence_curated_v6.xml','notebooks/Review_v6_EXECUTED.ipynb']:
 if not (R/f).is_file():missing.append(f)
assert not missing, 'Missing expected files: '+repr(missing[:30])
checks=R/'checksums.sha256'
if checks.exists():
 for line in checks.read_text().splitlines():
  sha, rel=line.split('  ',1)
  cur=hashlib.sha256((R/rel).read_bytes()).hexdigest()
  assert cur==sha, f'Checksum mismatch: {rel}'
print(f'v6 repository static checks passed: {len(rows)} supplemental IDs, 7 figures, critical models, checksums')
