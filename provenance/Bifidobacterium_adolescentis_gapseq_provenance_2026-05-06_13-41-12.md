# Bifidobacterium adolescentis GCF_000737885.1 gapseq reconstruction provenance

## Purpose

This document records the reconstruction and medium-specific gapfilling of a genome-scale metabolic model for *Bifidobacterium adolescentis* GCF_000737885.1. The model was reconstructed for pairwise community modeling with *Megasphaera* sp. MJR8396C under the revised anaerobic MJR8396C medium.

## Run metadata

- Run date: 2026-05-06_13-41-12
- Working directory: `/Users/wtscott/Documents/MATLAB/agora_megasphaera_workflow/MJR8396C_gapseq_package/Bifidobacterium_adolescentis_gapseq_run`
- Raw terminal log: `/Users/wtscott/Documents/MATLAB/agora_megasphaera_workflow/MJR8396C_gapseq_package/Bifidobacterium_adolescentis_gapseq_run/logs/Bifidobacterium_adolescentis_gapseq_run_2026-05-06_13-41-12.log`
- gapseq executable: `/Users/wtscott/gapseq/gapseq`
- Protein FASTA: `/Users/wtscott/Downloads/Bifidobacterium_adolescentis_genome_assembly/ncbi_dataset/data/GCF_000737885.1/protein.faa`
- Genome FASTA: `/Users/wtscott/Downloads/Bifidobacterium_adolescentis_genome_assembly/ncbi_dataset/data/GCF_000737885.1/GCF_000737885.1_ASM73788v1_genomic.fna`
- Gapfilling medium: `/Users/wtscott/Documents/MATLAB/agora_megasphaera_workflow/MJR8396C_gapseq_package/MJR8396C_candidate_anaerobic_medium_gapseq_revised.csv`
- Final RDS model: `Bifidobacterium_adolescentis_GCF_000737885_1_filled_MJR_anaerobic.RDS`
- Final SBML model: `Bifidobacterium_adolescentis_GCF_000737885_1_filled_MJR_anaerobic.xml`

## Commands used

```bash
cd "/Users/wtscott/Documents/MATLAB/agora_megasphaera_workflow/MJR8396C_gapseq_package/Bifidobacterium_adolescentis_gapseq_run"

"/Users/wtscott/gapseq/gapseq" find -p all -b 200 -m Bacteria "/Users/wtscott/Downloads/Bifidobacterium_adolescentis_genome_assembly/ncbi_dataset/data/GCF_000737885.1/protein.faa"

"/Users/wtscott/gapseq/gapseq" find-transport -b 200 "/Users/wtscott/Downloads/Bifidobacterium_adolescentis_genome_assembly/ncbi_dataset/data/GCF_000737885.1/protein.faa"

"/Users/wtscott/gapseq/gapseq" draft \
  -r protein-all-Reactions.tbl \
  -t protein-Transporter.tbl \
  -p protein-all-Pathways.tbl \
  -u 200 \
  -l 100

"/Users/wtscott/gapseq/gapseq" fill \
  -m protein-draft.RDS \
  -n "/Users/wtscott/Documents/MATLAB/agora_megasphaera_workflow/MJR8396C_gapseq_package/MJR8396C_candidate_anaerobic_medium_gapseq_revised.csv" \
  -b 100

mv protein.RDS "Bifidobacterium_adolescentis_GCF_000737885_1_filled_MJR_anaerobic.RDS"
mv protein.xml "Bifidobacterium_adolescentis_GCF_000737885_1_filled_MJR_anaerobic.xml"
```

## Key observed reconstruction output

The raw gapseq output reported:

- Input type: protein FASTA
- Predicted Gram stain: positive
- Number of pathways considered: 1929
- Number of reference sequences used for alignments: 154781
- ORF coverage: 25.04%
- Transporter prediction completed
- Gene-reaction list: 1091 unique genes
- Draft model construction: 1160 / 1160 reactions processed
- LP solver: glpk

## Medium-specific gapfilling result

The model was gapfilled on:

`MJR8396C_candidate_anaerobic_medium_gapseq_revised.csv`

Observed gapfilling summary:

### 1. Initial gapfilling

- Utilized candidate reactions: 130
- Added reactions: 125
- Added core reactions: 33
- Final growth rate: 0.3885231 h^-1

### 2. Biomass gapfilling using core reactions only

- Filled biomass components: 0
- Added reactions: 0
- Final growth rate: 0.3885231 h^-1

### 2b. Anaerobic biomass gapfilling using core reactions only

- Filled biomass components: 0
- Added reactions: 0
- Final growth rate: 0.3885231 h^-1

### 3. Energy-source gapfilling with core reactions only

- Filled components: 14
- Components: Glycerol-3-phosphate-e0, 2-Oxoglutarate-e0, GLCN-e0, D-Lactate-e0, Methanol-e0, CELB-e0, Melibiose-e0, Sorbitol-e0, Pyruvate-e0, Citrate-e0, Inosine-e0, NAD-e0, Xanthosine-e0, Xylan-b-1-4-e0
- Added reactions: 43
- Final growth rate: 0.4200339 h^-1

### 4. Potential metabolic-product gapfilling with core reactions only

- Filled components: 3
- Components: GLCN-e0, D-Ribose-e0, Cytidine-e0
- Added reactions: 6
- Final growth rate: 0.4200339 h^-1

## Uptake at limit

- Glycine: 2
- D-Glucose: 10
- L-Aspartate: 2
- L-Methionine: 2
- L-Cysteine: 2

## Top 10 predicted produced metabolites

Units reported by gapseq: mmol / (gDW * hr)

| Metabolite | Flux |
|---|---:|
| H+ | 20.27 |
| Acetate | 17.473 |
| CO2 | 7.19 |
| L-Lactate | 5.267 |
| H2O | 4.287 |
| NH3 | 2.568 |
| H2S | 1.967 |
| MTTL | 1.943 |
| Propionate | 1.775 |
| L-Lysine | 0.98 |

## Warnings and platform notes

The run was performed on macOS. The following messages were observed and should be reported as platform warnings rather than necessarily model-fatal errors:

- `Critical: (310.5) External MBEDTLS version mismatch: 3.6.5 headers vs. 3.6.6 runtime`
- `ps: illegal option -- q`
- `stat: illegal option -- c`

Despite these warnings, pathway prediction, transporter prediction, draft reconstruction, and gapfilling completed successfully.

## Interpretation for community modeling

The resulting *B. adolescentis* model is medium-compatible with the revised MJR8396C anaerobic medium after gapseq gapfilling, with final predicted growth of 0.4200339 h^-1. This is suitable as a co-culture-compatible reconstruction for MICOM simulations, but the model should be interpreted as a medium-gapfilled draft. Because 125 reactions were added during initial gapfilling and 43 additional reactions during energy-source gapfilling, downstream community claims should distinguish model-enabled growth from direct experimental validation.

