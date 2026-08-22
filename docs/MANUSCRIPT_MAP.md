# Manuscript-to-repository map

This map is aligned to the 21 August 2026 clean manuscript, **Genome-level constraint-based analysis of Megasphaera sp. MJR8396C reveals nutrient-dependent fermentation and lactate-linked cross-feeding**.

| Manuscript item | Repository source |
|---|---|
| Figure 1 — reconstruction and nutrient-dependent metabolic allocation | `figures/manuscript/Figure_1.{pdf,svg,png}`; `tables/Table_1_model_quality_and_audit_summary.csv`; `tables/Table_3_growth_and_exchange_three_media.csv` |
| Figure 2 — external physiological benchmarking | `figures/manuscript/Figure_2.{pdf,svg,png}`; `benchmarking/benchmark_model_predictions.csv`; `tables/Table_S5_external_benchmark_simulations.csv`; `tables/Table_S6_external_benchmark_interpretation.csv` |
| Figure 3 — comparative genomic context | `figures/manuscript/Figure_3.{pdf,svg,png}`; `tables/Table_S7_comparative_genomic_context.csv` |
| Figure 4 — substrate use, amino-acid interactions, butyrate feasibility | `figures/manuscript/Figure_4.{pdf,svg,png}`; `tables/Table_S11_carbon_normalized_substrate_screen.csv`; `tables/Table_S12_amino_acid_leave_one_out_revised.csv`; `tables/Table_S13_amino_acid_pair_simulations.csv`; `tables/Table_S14_butyrate_production_envelope.csv` |
| Figure 5 — quantitative StrainDesign predictions | `figures/manuscript/Figure_5.{pdf,svg,png}`; `tables/Table_4a_straindesign_solution_phenotypes.csv`; `tables/Table_4b_straindesign_intervention_sets.csv`; `tables/Table_S15_strain_design_phenotypes_with_wild_type.csv`; `tables/Table_S16_strain_design_target_frequency.csv`; `tables/Table_S17_wild_type_pfba_key_fluxes.csv`; `tables/Table_S18_Figure5_representative_designs.csv`; `tables/Table_S19_Figure5_relative_changes_percent.csv` |
| Figure 6 — pathway-centered StrainDesign interpretation | `figures/manuscript/Figure_6.{pdf,svg,png}` plus the StrainDesign tables above |
| Figure 7 — MICOM community growth and exchange | `figures/manuscript/Figure_7.{pdf,svg,png}`; `tables/Table_5_micom_exchange_fraction_0.95.csv`; `tables/Table_S20_micom_tradeoff_scan.csv`; `tables/Table_S21_micom_exchange_classification.csv`; `tables/Table_S22_micom_member_growth_fraction_0.95.csv`; `tables/Table_S23_micom_lactate_scfa_exchange_table.csv` |
| Table 1 — structural and QC properties | `tables/Table_1_model_quality_and_audit_summary.csv`; `analysis_summary.json` |
| Table 2 — gutSMASH pathway comparison | `tables/Table_2_gutsmash_pathway_comparison.csv`; `tables/Table_S8_gutSMASH_GEM_literature_triangulation.csv` |
| Table 3 — growth and key exchange fluxes | `tables/Table_3_growth_and_exchange_three_media.csv` |
| Table 4 — representative StrainDesign strategies | `tables/Table_4a_straindesign_solution_phenotypes.csv`; `tables/Table_4b_straindesign_intervention_sets.csv` |
| Table 5 — MICOM exchange fluxes at tradeoff 0.95 | `tables/Table_5_micom_exchange_fraction_0.95.csv` |
| Medium definitions | `media/`; `tables/Table_S1_medium_definitions_three_media.csv`; `tables/Table_S2_medium_to_exchange_mapping.csv` |
| Blocked reactions / non-GPR audit | `tables/Table_S3A_MEMOTE_blocked_reactions.tsv`; `tables/Table_S3B_gut-like-medium_blocked _reactions.tsv`; `tables/Table_S4_reactions_without_GPR.tsv`; `tables/Table_S4A_reactions_without_GPR_summary.tsv` |
| Model-quality reports | `quality_control/` |
| Reconstruction provenance | `provenance/` |

## Figure naming note

The uploaded project package contained several earlier figure-numbering variants. This GitHub-ready version places the **current manuscript-facing copies** in `figures/manuscript/` under canonical names `Figure_1` through `Figure_7` to avoid ambiguity. Numerical source tables are retained under their manuscript table identifiers.
