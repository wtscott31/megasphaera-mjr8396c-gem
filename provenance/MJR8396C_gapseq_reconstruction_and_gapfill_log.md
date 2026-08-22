# gapseq reconstruction and gap-filling provenance note

This file records the command-line reconstruction and gap-filling information supplied by the user for the latest Megasphaera sp. MJR8396C/MJR9396C gapseq rerun. It should be kept with the manuscript repository so that model provenance, medium assumptions, and gap-filling changes are auditable.

## Pathway/reaction search stage

```
Checking for pathways and reactions in:
/Users/wtscott/Documents/MATLAB/agora_megasphaera_workflow/MJR8396C_gapseq_package/2775506928_genes.faa

Pathways|Enzyme-Test|seed|kegg
Number of pathways to be considered: 1929

2 sequence files need to be downloaded from UniProt (via EC / metacyc-genes):
  1 / 2 (rev/1.2.1.59.fasta)
  2 / 2 (unrev/1.2.1.59.fasta)

Number of reference sequences used for alignments: 154781
ORF coverage: 25.81 %
Running time: 160 s.
Protein fasta detected.
```

The run emitted external library/runtime and macOS compatibility messages:

```
Critical: (310.5) External MBEDTLS version mismatch: 3.6.5 headers vs. 3.6.6 runtime
ps: illegal option -- q
stat: illegal option -- c
```

These messages should be documented as runtime warnings. They do not by themselves prove model failure, but they matter for reproducibility because they indicate non-Linux command assumptions inside the workflow and an MBEDTLS runtime/header mismatch.

## Transporter stage

The transporter screen detected transport capacity for organic acids, sugars, amino acids, ions, cofactors, bile-related compounds, polyamines, nucleosides, and stress-related metals. Notable detected substrates include acetate, propionate, D-/L-lactate, formate, fumarate, succinate, pyruvate, glucose, xylose, amino acids, GABA, taurine, thiosulfate, sulfate, thiamin, riboflavin, pantothenate, biotin, vitamin B12, and multiple metal ions.

Transporters without corresponding gapseq database reactions were reported for cholesterol, cyanate, cytidine, FAD, L-fucose, L-rhamnose, nitrite, turanose, xanthosine, and xylan-b-1-4. Alternative transport reactions were added for cholesterol, cyanate, cytidine, FAD, L-arabinol, L-fucose, L-rhamnose, pyridoxal, turanose, and xanthosine.

## Draft model construction

```
Creating Gene-Reaction list... 1357 unique genes
Constructing draft model:
 1428 / 1428
LP solver: glpk
Loading model file 2775506928_genes-draft.RDS
```

## Anaerobic medium gap-fill

Command:

```
~/gapseq/gapseq fill \
  -m 2775506928_genes-draft.RDS \
  -n /Users/wtscott/Documents/MATLAB/agora_megasphaera_workflow/MJR8396C_gapseq_package/MJR8396C_candidate_anaerobic_medium_gapseq_revised.csv
```

Gap-fill output:

```
1. Initial gapfilling:
Utilized candidate reactions: 112
Added reactions: 107
Added core reactions: 64
Final growth rate: 0.6939079

2. Biomass gapfilling using core reactions only:
Filled components: 2 (TPP-c0, CoA-c0)
Added reactions: 14
Final growth rate: 0.6939079

2b. Anaerobic biomass gapfilling:
Filled components: 0
Added reactions: 14
Final growth rate: 0.6939079

3. Energy-source gapfilling:
Filled components: 17
(L-Lysine-e0, Glycerol-3-phosphate-e0, L-Valine-e0, Guanosine-e0, L-Isoleucine-e0, L-Tryptophan-e0, Xylose-e0, L-Idonate-e0, Methanol-e0, Inosine-e0, D-Galacturonate-e0, D-Ribose-e0, NAD-e0, 2-Oxoglutarate-e0, Adenosine-e0, Cytidine-e0, Xanthosine-e0)
Added reactions: 70
Final growth rate: 0.7000702

4. Product-potential gapfilling:
Filled components: 7
(Ornithine-e0, GLCN-e0, D-Lactate-e0, Taurine-e0, Neu5Ac-e0, D-Ribose-e0, Cytosine-e0)
Added reactions: 29
Final growth rate: 0.7004949

Uptake at limit:
D-Glucose:10, L-Aspartate:2, L-Methionine:2, L-Cysteine:2, Acetate:5, Glycine:2, L-Lactate:10

Top 10 produced metabolites [mmol/(gDW*h)]:
H2O:24.854, CO2:21.733, Propionate:11.261, Butyrate:10.847, H+:4.619, H2S:1.932, MTTL:1.893, NH3:1.472, Succinate:0.012, 5-Methylthio-D-ribose:0.002
```

## Vaginal/BV-like medium gap-fill

Command:

```
~/gapseq/gapseq fill \
  -m 2775506928_genes-draft.RDS \
  -n /Users/wtscott/Documents/MATLAB/agora_megasphaera_workflow/MJR8396C_gapseq_package/MJR8396C_candidate_vaginal_BV_like_medium_gapseq_revised.csv
```

Gap-fill output:

```
1. Initial gapfilling:
Utilized candidate reactions: 114
Added reactions: 108
Added core reactions: 65
Final growth rate: 0.4698331

2. Biomass gapfilling using core reactions only:
Filled components: 2 (TPP-c0, L-Glutamate-c0)
Added reactions: 15
Final growth rate: 0.4698331

2b. Anaerobic biomass gapfilling:
Filled components: 0
Added reactions: 15
Final growth rate: 0.4698331

3. Energy-source gapfilling:
Filled components: 15
(L-Lysine-e0, Glycerol-3-phosphate-e0, L-Valine-e0, L-Isoleucine-e0, Xylose-e0, L-Idonate-e0, Methanol-e0, Inosine-e0, D-Galacturonate-e0, D-Ribose-e0, NAD-e0, 2-Oxoglutarate-e0, Adenosine-e0, Cytidine-e0, Xanthosine-e0)
Added reactions: 82
Final growth rate: 0.4743943

4. Product-potential gapfilling:
Filled components: 8
(Ornithine-e0, Spermidine-e0, GLCN-e0, Taurine-e0, Neu5Ac-e0, D-Ribose-e0, Cytosine-e0, PAN-e0)
Added reactions: 21
Final growth rate: 0.5468209

Uptake at limit:
L-Glutamate:5, D-Glucose:5, L-Methionine:1, L-Cysteine:1, L-Serine:5, Pyruvate:2, Acetate:2, D-Lactate:10

Top 10 produced metabolites [mmol/(gDW*h)]:
Formate:20.477, H+:16.088, Butyrate:11.446, CO2:7.263, NH3:6.397, H2O:4.694, Succinate:4.416, MTTL:1.204, Propionate:1.004, H2S:0.659
```

## Recommended manuscript wording

The rerun should be described as a gapseq v2-series reconstruction/gap-fill based on protein FASTA `2775506928_genes.faa` and the two machine-readable media files. Because the SBML model notes report `gapseq version: 2.0.0 c64db99b`, while the rerun was described as using the v2.0.1 release, the final manuscript and repository should explicitly reconcile the exact executable commit or release used to create the deposited SBML.
