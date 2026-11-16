# Reproduction and software environment

## Starting point

The original model is `models/original/Megasphaera_sp_MJR8396C_filled_anaerobic.xml` and the original partner model is `models/original/Bifidobacterium_adolescentis_filled_anaerobic.xml`. The exact original reconstructions, sequence identifiers, input media, GPR adjustments and all accepted v6 summary outputs are deposited separately. Never silently overwrite the original with a curated variant.

## Environment

The execution used Python 3.11 with COBRApy, MICOM 0.39.1, highspy 1.15.1, Biopython, and BOFdat pinned to commit `cf2de074c2d789ad47b890083cbbffadab40f177`; see `requirements-resolved.txt`. Native uniform matched reconstruction used gapseq 2.1.0 and Bacteria DB version 1.5 (Zenodo record 20446806). gapseq, BOFdat, and other third-party code and large reference database files must be obtained from their original vendors.

To install a Python environment:

```bash
python3.11 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-resolved.txt
python scripts/verify_repository_v6.py
```

The pinned resolved requirements include optional commercial-solver packages and a Git dependency; adapt installation to your platform and licensing, recording all changes. The native HiGHS / custom QP implementation is explicitly differentiated from a stock commercial-solver run.

## Computational analysis

The original execution sequence is preserved in `docs/ORIGINAL_EXECUTION_ORDER.txt`. The matched gapseq workflow is in `scripts/reconstruct_matched_v6.py`; set `GAPSEQ_BIN` and `PATH` as needed. The workflow expects R and GLPK and the Bacteria 1.5 database. Rerun of the full matched MICOM grids requires corresponding full input and output files from the Zenodo archive. `notebooks/Review_v6_EXECUTED.ipynb` is an executed review of those files. Its outputs are retained, but this does not mean every full-grid test executes in isolation from the external database/cache.

## Provenance / outputs

`Supplementary_Data S1–S87` map original data and accepted revision analyses into canonical release paths. The `supplementary/SUPPLEMENTARY_INDEX.csv` is the authoritative file map. Key activities included reciprocal ANI and GTDB reference lookup; only the exact assembly GTDB lookup is claimed, not a newly executed phylogenetic GTDB-Tk placement. Biomass predictions are sensitivity calculations and do not represent measured target-strain biomass fractions.

The source dataset also includes curated-comparator scenarios, not fully experimentally calibrated cross-species predictions. No direct strain-specific growth or exchange measurements were created.
