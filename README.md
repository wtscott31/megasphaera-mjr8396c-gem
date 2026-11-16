# Genome-scale modeling of *Megasphaera* sp. MJR8396C

**Companion code and source data for:** “Genome-scale constraint-based modeling of *Megasphaera* sp. MJR8396C predicts nutrient-dependent fermentation and conditional microbial exchange.”

This repository contains the authors' computational revision based on the deposited MJR8396C model and uniform gapseq 2.1.0 reconstructions of MJR8396C, *Megasphaera indica* NM10, *Megasphaera* BL7 and *Bifidobacterium adolescentis* 22L. It is **not** a claim of direct experimental validation of MJR8396C, nor does it establish intestinal colonization: the reference isolate derives from a human vaginal sample.

The repo provides immutable original SBML files, derived sequence-supported GPR models, machine-readable source tables S1–S87, accession-resolved genome FASTA, executed notebook, reproducible code and original input provenance. Larger full-grid community flux matrices, complete matched reconstruction intermediate files, and diagnostic model variants are prepared separately for a Zenodo data archive. **No Zenodo DOI has been minted** by preparing this package.

## Main findings and interpretive constraints

- The original 1,849-reaction, 1,614-metabolite, 477-gene model is a draft, gap-filled reconstruction; GEM/MEMOTE consistency is distinct from physiological calibration.
- Growth and secretion depend on assumed medium, biomass and maintenance constraints; pFBA points are not unique or obligatory flux predictions.
- Targeted Bcd–Etf, Rnf and PFOR gene associations are supported by sequence comparisons; other electron-transfer and terminal-pathway assignments remain uncertain.
- Historical reaction-only nucleotide-associated butyrate interventions are not feasible gene deletions under the deposited full-network GPR logic.
- Modeled lactate exchange is conditional. Blocking bifidobacterial lactate export does not demonstrate an obligate lactate dependency. Equal-growth MICOM outcomes can reflect its cooperative-tradeoff regularization.

## Navigate

| Folder | Contents |
|---|---|
| `models/original` | Original deposited SBML models (unchanged) |
| `models/` | Sequence-annotated matched reconstructions and named sensitivity models |
| `inputs/genomes` | Accession-verified genome FASTA, proteins and assembly provenance |
| `inputs/submission_reference` | Original manuscript medium definitions and source tables needed by legacy workflows (not previous version releases) |
| `scripts` / `tests` | Modeling, quality checks, and numerical solver adapters |
| `notebooks` | Executed v6 results overview |
| `supplementary/SUPPLEMENTARY_INDEX.csv` | Canonical S1–S87 source-data map with exact paths |
| `figures/manuscript` | Seven revised manuscript figures and supplementary sequence map |
| `results/current_summary` | Current numerical results and test certificates |
| `provenance` | Versioned methods, input metadata, and source hashes |

## Reproducibility

See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) and [`docs/RESULTS_AND_LIMITATIONS.md`](docs/RESULTS_AND_LIMITATIONS.md). Run `python scripts/verify_repository_v6.py` to verify presence, supplementary paths and SHA-256 checksums. Full reconstruction requires separately installed pinned gapseq with database, R/GLPK and the full Zenodo intermediate data; the historical environment used a Mac-specific Conda location, which is not bundled. For a complete rerun, use the recorded versioned external tools and command sequence rather than treating a notebook as a one-command reproduction.

## Data deposition and citation

Proposed GitHub URL: https://github.com/wtscott31/megasphaera-mjr8396c-gem (please push the versioned release after author approval). Intended Zenodo archive name: `megasphaera-mjr8396c-gem_v6_full_data.zip`. Do not insert a Zenodo DOI into the manuscript until an actual Zenodo record exists. The `CITATION.cff` is a software/data citation template, not an accepted article citation.

## Licensing

<<<<<<< Updated upstream
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

MIT License
=======
**No new license is assigned.** The authors should select and document an approved software license before making the code public. External database/reference models retain their respective terms and are attributed in the provenance tables. Vendored upstream gapseq and BOFdat implementations are intentionally excluded.
>>>>>>> Stashed changes
