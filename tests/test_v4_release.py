from pathlib import Path
import json,sys
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))

def test_all_deposited_genes_map_by_exact_sequence_or_coordinates():
    d=pd.read_csv(ROOT/'results/JGI_mapping/deposited_model_477_gene_mapping.csv')
    assert len(d)==477 and d.JGI_gene.nunique()==477
    assert d.public_protein_accessions.notna().all()
    assert d.public_locus_tags.notna().all()
    assert d.original_sequence_verified_in_original_genome.all()
    assert d.original_sequence_verified_in_public_genome.all()
    assert set(d.mapping_status)=={'unique_exact_protein_match','unique_coordinate_verified_annotation_boundary_variant'}

def test_original_contigs_are_exact_public_segments_but_not_same_assembly():
    d=pd.read_csv(ROOT/'results/JGI_mapping/exact_contig_comparison.csv');s=json.loads((ROOT/'results/JGI_mapping/summary.json').read_text())
    assert len(d)==80 and d.exact_match.all()
    assert s['original_total_bp']==2426262 and s['public_total_bp']==2453897
    assert not s['identical_assembly']

def test_annotation_only_model_preserves_all_chemistry_bounds_and_GPR():
    from native_analysis import load,ROOT as root
    a=load();b=load(path=root/'models/MJR8396C_sequence_annotated_primary.xml')
    for r in a.reactions:
        q=b.reactions.get_by_id(r.id)
        assert r.bounds==q.bounds and r.gene_reaction_rule==q.gene_reaction_rule
        assert {m.id:v for m,v in r.metabolites.items()}=={m.id:v for m,v in q.metabolites.items()}
    assert all(g.annotation.get('ncbiprotein') for g in b.genes)

def test_original_JGI_native_BOFdat_failure_and_transparent_correction():
    d=json.loads((ROOT/'results/JGI_mapping/annotation_and_original_DNA_audit.json').read_text())
    assert abs(d['native_fraction_recovered']-.5)<1e-12
    assert abs(sum(d['genome_derived_double_strand_frequencies'].values())-1)<1e-12
    assert abs(d['baseline_net_mass']-d['candidate_net_mass'])<1e-10
    assert d['native_raw_rejected']

def test_related_strain_empirical_candidate_is_mass_closed_and_honest():
    d=json.loads((ROOT/'results/empirical_biomass/candidate_assumptions.json').read_text())
    assert d['candidate_only'] and not d['target_strain_measurements_available']
    assert abs(d['before_total_mass']-d['after_total_mass'])<1e-10
    p=pd.read_csv(ROOT/'results/empirical_biomass/Wallace1986_measured_amino_acid_pools.csv')
    assert len(p)==16 and p.determinations.eq(3).all()

def test_empirical_candidate_introduces_no_new_precursor_failure():
    p=pd.read_csv(ROOT/'results/empirical_biomass/precursor_demands.csv')
    blocked={g:set(d.loc[d.demand_max<1e-6,'metabolite']) for g,d in p.groupby('model')}
    assert len(blocked)==3 and len({frozenset(s) for s in blocked.values()})==1

def test_empirical_candidate_does_not_create_energy():
    d=pd.read_csv(ROOT/'results/empirical_biomass/ATP_cycle_checks.csv')
    assert len(d)==6 and d.status.eq('optimal').all() and d.ATP_max.abs().max()<1e-8

def test_all_240_new_native_MICOM_cases_pass_full_network_checks():
    d=pd.read_csv(ROOT/'results/native_micom_v4/native_tradeoff_controls.csv')
    assert len(d)==240 and d.status.eq('optimal').all()
    for col in ['mass_residual','bound_violation','QP_feasibility','QP_stationarity','QP_complementarity']:
        assert d[col].max()<1e-6
    assert d.QP_global_gap.max()<1e-7

def test_new_native_MICOM_matches_independent_v2():
    d=pd.read_csv(ROOT/'results/native_micom_v4/native_vs_independent_v2.csv')
    assert d.growth_M_difference.abs().max()<2e-6
    assert d.growth_B_difference.abs().max()<2e-6
    assert d.community_max_difference.abs().max()<1e-8

def test_new_native_exchange_accounting_and_donor_controls():
    d=pd.read_csv(ROOT/'results/native_micom_v4/native_taxon_community_fluxes.csv')
    assert d.accounting_residual.abs().max()<1e-6
    x=d[d.control.isin(['block_both','no_lactate_block_both'])&d.compound.isin(['cpd00159','cpd00221'])]
    assert x.B_specific.max()<1e-7

def test_QP_backend_analytic_and_full_rank_independent_cases():
    from conditional_gradient_qp import install
    from optlang.hybrid_interface import Model,Variable,Objective,Constraint
    install()
    cases=[(np.array([[2.,1.],[1.,4.]]),np.array([.8,1.2])),(np.diag([2.,4.,6.]),np.array([1.,2.,3.]))]
    for H,target in cases:
        m=Model();xs=[Variable('x'+str(i),lb=0,ub=10) for i in range(len(target))]
        m.add(xs+[Constraint(sum(xs),lb=float(target.sum()),ub=float(target.sum()))])
        expr=sum(.5*H[i,j]*xs[i]*xs[j] for i in range(len(xs)) for j in range(len(xs)))-sum(float((H@target)[i])*xs[i] for i in range(len(xs)))
        m.objective=Objective(expr,direction='min');assert m.optimize()=='optimal'
        assert np.max(np.abs(np.array([v.primal for v in xs])-target))<1e-5
        assert m.problem.qp_certificates[-1]['global_suboptimality_bound']<1e-7

def test_new_native_QP_conditioned_FVA_has_finite_ordered_ranges():
    d=pd.read_csv(ROOT/'results/native_micom_v4/native_QP_conditioned_exchange_FVA.csv')
    assert len(d)==540 and not d[['minimum','maximum']].isna().any().any()
    assert (d.minimum<=d.maximum+1e-6).all()

def test_comparator_native_files_not_silently_repaired_to_force_agreement():
    d=json.loads((ROOT/'results/comparators/v4_iME375_supplement_discrepancies.json').read_text())
    assert d['actual_biomass_synthesis']=='R615' and d['deposited_R610_serine']==1.141
    assert d['S1_serine_mmol_per_g_protein']==.357
