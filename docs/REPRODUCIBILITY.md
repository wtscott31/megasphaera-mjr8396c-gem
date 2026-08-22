# Reproducibility notes

## What can be rerun directly

`notebooks/MJR8396C_revised_manuscript_simulations.ipynb` and `scripts/run_core_simulations.py` contain executable code for the core single-strain constraint-based analyses, including medium application, FBA/pFBA, the three-medium comparison, carbon-normalized substrate screening, amino-acid perturbations, and the butyrate production envelope.

The standalone script intentionally parses the objective from the deposited gapseq SBML rather than assuming a hard-coded biomass reaction name.

## What is deposited as finalized analysis output

The supplied package includes finalized StrainDesign and MICOM result tables and manuscript figures. The general notebook contains an optional StrainDesign section and code that integrates existing MICOM outputs, but it does **not** contain a complete one-command rerun of the full dedicated MICOM workflow. Consequently, this repository should not claim a fully containerized, end-to-end regeneration of every manuscript panel from raw genome sequence.

This distinction is important for transparent reproducibility: the final numerical outputs and models are deposited, while the complete specialized execution history is only partially represented in the supplied scripts.

## Reconstruction provenance

The manuscript-facing provenance identifies gapseq **v2.0.0, commit `c64db99b`**, with bacterial sequence database v1.4. The provenance files are retained verbatim in `provenance/`.

## Model limitations relevant to reuse

The MJR8396C reconstruction is medium-gapfilled and incompletely curated. The manuscript reports 766 condition-specific blocked reactions in the gut-like medium, 677 reactions blocked under MEMOTE's standardized topology configuration, and 514 reactions without GPR support. Users should therefore avoid treating reaction-level engineering predictions or absolute extracellular fluxes as experimentally validated quantities.

## Recommended release practice

For a formal archival release, create a GitHub release matching the manuscript version and archive that release with Zenodo (or an equivalent repository) to obtain a DOI. Replace the manuscript's `[FINAL REPOSITORY RELEASE URL/DOI]` placeholder only after the immutable release/DOI exists.
