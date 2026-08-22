# Megasphaera sp. MJR8396C genome-scale metabolic modeling

[![Reproducibility check](https://github.com/USERNAME/REPOSITORY/actions/workflows/verify.yml/badge.svg)](https://github.com/USERNAME/REPOSITORY/actions/workflows/verify.yml)

Computational models, media definitions, quality-control reports, analysis notebooks, machine-readable source data, and manuscript figures supporting:

> **Genome-level constraint-based analysis of Megasphaera sp. MJR8396C reveals nutrient-dependent fermentation and lactate-linked cross-feeding**  
> Joy Iheoma van Geerestein*, Andrea Dell'Olio*, Josep Rubert, William T. Scott Jr.  
> *Co-first authors

The study reconstructs a medium-conditioned genome-scale metabolic model (GEM) of **Megasphaera sp. MJR8396C**, evaluates nutrient-dependent fermentation, benchmarks selected model predictions against published *Megasphaera* physiology, examines genomic pathway context, identifies candidate intervention strategies with StrainDesign, and models a two-member interaction with *Bifidobacterium adolescentis* using MICOM.

## Key model and study identifiers

- **MJR8396C genome assembly:** `GCA_001546855.1`
- **B. adolescentis community-model genome:** `GCF_000737885.1`
- **MJR8396C model:** `models/Megasphaera_sp_MJR8396C_filled_anaerobic.xml`
- **B. adolescentis model:** `models/Bifidobacterium_adolescentis_filled_anaerobic.xml`
- **gapseq executable provenance:** v2.0.0, commit `c64db99b`; bacterial sequence database v1.4
- **Final MJR8396C GEM:** 1,849 reactions, 1,614 metabolites, 477 genes
- **MEMOTE total score:** 88.2%

The model is a **medium-gapfilled reconstruction intended for hypothesis generation**, not a quantitatively calibrated digital twin. The manuscript explicitly reports remaining uncertainty, including blocked reactions, non-GPR reactions, and incomplete quantitative agreement with lactate physiology.

## Repository structure

```text
.
├── models/              # Final SBML models
├── media/               # Gut-like, carbohydrate-rich, and amino-acid-rich media
├── notebooks/           # Primary and executed analysis notebooks
├── scripts/             # Standalone core simulation script + verification utility
├── tables/              # Machine-readable data behind main/supplementary tables
├── figures/
│   ├── manuscript/      # Canonical current Figures 1-7 (PDF/SVG/PNG)
│   └── supplementary/   # Supplementary figure(s)
├── benchmarking/        # External physiological benchmarking outputs
├── quality_control/     # MEMOTE reports and model-comparison reports
├── provenance/          # gapseq reconstruction/gapfilling logs
├── docs/                # Manuscript-to-repository and reproducibility notes
└── .github/workflows/   # Lightweight repository integrity check
```

## Quick start

A lightweight standalone implementation reproduces the core FBA/pFBA analyses using SciPy/HiGHS and reads the objective directly from the deposited SBML.

```bash
python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python scripts/run_core_simulations.py
```

For notebook-based analyses:

```bash
jupyter lab notebooks/MJR8396C_revised_manuscript_simulations.ipynb
```

The notebook uses COBRApy for model analysis and contains an **optional** StrainDesign section. A MILP-capable solver is required for StrainDesign. The deposited MICOM and StrainDesign result tables/figures are included in the repository; the current general notebook integrates those outputs rather than rerunning the full dedicated MICOM workflow.

## Manuscript-facing analyses

The deposited package covers the manuscript's computational workflow:

- model reconstruction and quality audit;
- three-medium pFBA/FVA analysis;
- external physiological benchmarking;
- gutSMASH/GEM comparison source tables;
- carbon-normalized substrate screening;
- amino-acid leave-one-out and pair-interaction analyses;
- butyrate production envelope;
- StrainDesign intervention results;
- MICOM cooperative-tradeoff and taxon-resolved exchange results.

See [`docs/MANUSCRIPT_MAP.md`](docs/MANUSCRIPT_MAP.md) for the exact mapping from manuscript figures/tables to repository files.

## Reproducibility scope

The repository contains the final SBML models, machine-readable media, model-quality reports, numerical result tables, current manuscript figures, and the core executable notebook/script. Some specialized analyses (particularly the full MICOM and StrainDesign runs) are represented by their deposited outputs and integration code rather than a single end-to-end rerun script. This is stated explicitly in [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).

To run a fast integrity check:

```bash
python scripts/verify_repository.py
```

## Interpretation guardrails

1. Predicted biomass fluxes are model outputs under specified exchange bounds and should not be described as measured growth rates.
2. Positive exchange flux denotes secretion and negative flux denotes uptake.
3. MICOM results are condition-specific, model-predicted cross-feeding hypotheses.
4. StrainDesign reaction deletions are network-level intervention hypotheses, not validated gene-editing designs.
5. The current reconstruction should not be treated as quantitatively validated for strain-specific extracellular metabolomics.

## Citation

Please cite the associated manuscript/preprint and the software used in the workflow (gapseq, MEMOTE, COBRApy, StrainDesign, and MICOM) as appropriate. A machine-readable citation template is provided in [`CITATION.cff`](CITATION.cff).

## License

No repository license has been assigned in the supplied project files. **Before making the repository public, select and add an appropriate license for the code, models, and data.** This is intentionally left unresolved rather than assigning legal terms without author approval.
