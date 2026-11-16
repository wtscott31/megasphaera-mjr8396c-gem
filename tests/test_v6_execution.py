from pathlib import Path
import json
import numpy as np
import pandas as pd
import pytest
import cobra

ROOT=Path(__file__).resolve().parents[1]
cobra.Configuration().solver='glpk'

@pytest.mark.parametrize('strain',['MJR8396C','NM10','BL7','B22L'])
def test_matched_pipeline_completed_with_valid_models(strain):
    d=ROOT/'results/matched_reconstruction_v6'/strain
    for stage in ['find','transport','draft','fill']:
        s=json.loads((d/(stage+'_status.json')).read_text());assert s['returncode']==0 and s['expected_output_exists']
    m=cobra.io.read_sbml_model(str(d/(strain+'.xml')))
    assert len(m.reactions)>1000 and len(m.genes)>300
    assert 'Sequence DB version: 1.5' in (d/(strain+'-all-Reactions.tbl')).read_text().splitlines()[1]

@pytest.mark.parametrize('pair',['MJR8396C','NM10','BL7','GPRcurated_original22L','GPRcurated_filtered22L','MJR8396C_matched_curated_sensitivity','NM10_matched_curated_sensitivity','BL7_matched_curated_sensitivity'])
def test_complete_native_community_and_certificates(pair):
    d=ROOT/'results/matched_micom_v6'/pair;f=pd.read_csv(d/'native_tradeoff_controls.csv')
    assert len(f)==240 and set(f.status)=={'optimal'}
    assert f.mass_residual.max()<1e-6 and f.bound_violation.max()<1e-6
    assert f.QP_stationarity.max()<1e-6 and f.QP_global_gap.max()<1e-6
    assert np.max(abs(f.community_growth-f.weighted_growth))<1e-6
    p=pd.read_csv(d/'native_taxon_community_fluxes.csv');assert p.accounting_residual.abs().max()<1e-6
    assert np.isfinite(p[['M_specific','B_specific','community']].values).all()
    fv=pd.read_csv(d/'exchange_FVA95.csv');assert len(fv)==540 and 'error' not in fv
    assert (fv.minimum<=fv.maximum+1e-6).all()

@pytest.mark.parametrize('strain',['MJR8396C','NM10','BL7','B22L'])
def test_full_AA_maintenance_and_loopless_checks(strain):
    d=ROOT/'results/matched_monocultures_v6'/strain
    aa=pd.read_csv(d/'AA_singles_and_pairs.csv');assert len(aa)==210
    assert aa.AA2.notna().sum()==190
    assert len(pd.read_csv(d/'maintenance_uptake_grid.csv'))==108
    f=pd.read_csv(d/'loopless_FVA95_hybrid.csv')
    assert set(f.medium)=={'gut_like','carbon_rich','AA_rich'} and set(f.status)=={'optimal'}
    assert np.isfinite(f[['minimum','maximum']].values).all()
    assert (f.minimum<=f.maximum+1e-6).all()
    en=pd.read_csv(d/'ATP_cycle_checks.csv');assert en.ATP_max.abs().max()<1e-6

def test_primary_GPR_curation_preserves_every_numeric_constraint():
    a=cobra.io.read_sbml_model(str(ROOT/'models/MJR8396C_sequence_annotated_primary.xml'))
    b=cobra.io.read_sbml_model(str(ROOT/'models/MJR8396C_sequence_supported_GPR_curated_v6.xml'))
    changes=[]
    for old,new in zip(a.reactions,b.reactions):
        assert old.id==new.id and old.bounds==new.bounds
        assert {m.id:v for m,v in old.metabolites.items()}=={m.id:v for m,v in new.metabolites.items()}
        if old.gene_reaction_rule!=new.gene_reaction_rule:changes.append(old.id)
    assert set(changes)=={'rxn90001_c0','rxn90053_c0','rxn05938_c0'}
    assert len(b.genes)==483
    assert len(b.reactions.rxn90053_c0.genes)==6 and len(b.reactions.rxn90001_c0.genes)==3

def test_22L_evidence_and_conservative_model_scope():
    f=pd.read_csv(ROOT/'results/functional_curation_v6/B22L_deposited_gene_accession_crosscheck.csv')
    assert len(f)==349 and f.present_exact_accession_version.all()
    s=json.loads((ROOT/'results/functional_curation_v6/B22L_inventory.json').read_text())
    x={r['component']:r for r in s};assert x['L_NAD_LDH']['candidate_proteins']=='WP_003810784.1'
    assert x['D_NAD_LDH']['candidate_proteins']=='' and x['LarA']['candidate_proteins']==''
    m=cobra.io.read_sbml_model(str(ROOT/'models/B22L_D_LDH_unsupported_gapfill_disabled_v6.xml'))
    assert m.reactions.rxn00500_c0.bounds==(0,0)
    f=pd.read_csv(ROOT/'results/matched_micom_v6/GPRcurated_filtered22L/native_taxon_community_fluxes.csv')
    assert f[f.compound.eq('cpd00221')].B_specific.abs().max()<1e-6

def test_NM10_frameshift_is_traceable_not_assumed_absence():
    s=json.loads((ROOT/'results/matched_functional_curation_v6/NM10_Bcd_all_feature_check.json').read_text())
    assert any('potential frameshift' in str(r['qualifiers']) for r in s)
    f=pd.read_csv(ROOT/'results/matched_functional_curation_v6/NM10_Bcd_high_identity_fragments.csv')
    assert len(f)>=2 and f.pident.min()>90 and f.query_coverage_percent.max()<70
    e=pd.read_csv(ROOT/'results/matched_functional_curation_v6/matched_GPR_and_bound_edits.csv')
    assert 'STRICT SENSITIVITY' in e[(e.model=='NM10')&(e.reaction=='rxn90001_c0')].iloc[0].reason

def test_intervention_gene_counterfactuals_remain_nonviable():
    f=pd.read_csv(ROOT/'tables/S78_matched_interventions_v6.csv')
    x=f[f.model.isin(['MJR8396C','NM10','BL7'])&f.design.isin(['OK2','OK3','OK4'])&f.level.eq('gene')]
    assert len(x)==18 and x.growth.abs().max()<1e-7

def test_common_carbon_budget_and_cycle_correction():
    d=ROOT/'results/comparator_calibration_v6';f=pd.read_csv(d/'common_carbon_budget_simulations.csv')
    x=f[f.status.eq('optimal')];assert np.isfinite(x.actual_gross_organic_C_uptake).all()
    assert x.actual_gross_organic_C_uptake.max()<=60+1e-6
    en=pd.read_csv(d/'comparator_ATP_cycle_checks.csv')
    assert en[en.model.eq('iME375')].max_ATP_drain.min()>100
    assert en[en.model.isin(['iME375_cycle_constrained','iME375_partial_plus_cycle'])].max_ATP_drain.abs().max()<1e-6
    norm=pd.read_csv(d/'formula_mass_normalization_decisions.csv')
    assert norm.status.eq('formula_normalized_candidate').sum()==4
    assert norm.status.eq('excluded_from_mass_normalization').sum()==2

def test_biomass_candidates_preserve_mass_and_make_precursors():
    d=ROOT/'results/matched_biomass_v6';f=pd.read_csv(d/'growth_maintenance_sensitivity.csv')
    assert len(f)==216 and f.biomass_mass.between(.9999,1.0001).all()
    assert pd.read_csv(d/'ATP_cycle_checks.csv').ATP_max.abs().max()<1e-6
    p=pd.read_csv(d/'precursor_producibility.csv')
    assert set(p[p.max_demand.lt(1e-7)].metabolite)=={'cpd11493_c0','cpd00166_c0','cpd15665_c0'}
    paired=pd.read_csv(d/'paired_carrier_precursor_demands.csv');assert len(paired)==27 and paired.max_paired_demand.min()>1e-7
    scope=json.loads((ROOT/'results/comparator_calibration_v6/iME375_partial_repair_scope.json').read_text())
    assert scope['whole_cell_mass_resolved'] is False

def test_effective_supply_comparison_reports_unattainable_growth_honestly():
    x=pd.read_csv(ROOT/'tables/S80b_effective_supply_mono_community_v6.csv');assert len(x)==120
    missing=x.matched_growth_FVA_status.eq('unmatched_growth')
    assert (x.loc[missing,'community_member_growth']>x.loc[missing,'monoculture_growth_max']+1e-6).all()
    valid=x[~missing]
    assert np.isfinite(valid[['mono_butyrate_min_at_matched_member_growth','mono_butyrate_max_at_matched_member_growth']].values).all()
    assert (valid.mono_butyrate_min_at_matched_member_growth<=valid.mono_butyrate_max_at_matched_member_growth+1e-6).all()
