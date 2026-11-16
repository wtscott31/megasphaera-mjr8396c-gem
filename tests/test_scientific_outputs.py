from pathlib import Path
import hashlib,json,zipfile,sys
import numpy as np,pandas as pd,pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))

def test_original_source_models_match_input_archive():
    with zipfile.ZipFile(ROOT/'provenance/MJR8396C_executed_revision_v2_COMPLETE.zip') as z:
        for p in (ROOT/'baseline_v2/megasphaera-mjr8396c-gem_major-revision/models').glob('*.xml'):
            names=[n for n in z.namelist() if n.endswith('/models/'+p.name) and '/extended/' not in n]
            assert len(names)==1
            assert hashlib.sha256(z.read(names[0])).digest()==hashlib.sha256(p.read_bytes()).digest()
def test_downloaded_sequences_pass_NCBI_md5():
    f=pd.read_csv(ROOT/'results/NCBI_md5_validation.csv');assert len(f)>10;assert f.matches.all()
def test_genome_statistics_against_independent_NCBI_metadata():
    f=pd.read_csv(ROOT/'results/genome_statistics.csv').set_index('accession');r=f.loc['GCA_001546855.1']
    meta=json.loads((ROOT/'provenance/assemblies/GCA_001546855.1.json').read_text())['reports'][0]
    assert r.total_bp==int(meta['assembly_stats']['total_sequence_length'])
    assert r.canonical_bp==int(meta['assembly_stats']['atgc_count'])
    assert r.total_bp-r.canonical_bp==r.ambiguous_bp==5500
    assert r.scaffold_N50_bp==meta['assembly_stats']['scaffold_n50']
    assert r.N_split_contig_N50_bp==meta['assembly_stats']['contig_n50']
    assert r.features_CDS==r.protein_fasta_records==2422
    assert abs(r.GC_percent_canonical-100*(r['C']+r['G'])/(r['A']+r['C']+r['G']+r['T']))<1e-10
def test_GenBank_RefSeq_sequence_identity():
    assert json.loads((ROOT/'results/target_assembly_validation.json').read_text())['genbank_refseq_same_sequence_multiset']
def test_taxonomy_method_is_exact_lookup():
    t=json.loads((ROOT/'results/GTDB_exact_assembly_lookup.json').read_text())
    assert 'NOT a new GTDB-Tk' in t['method'];assert t['GTDB_Tk_status']=='not_run'
    assert t['result'][0][0]=='RS_GCF_001546855.1'
def test_no_absent_ANI_reported_as_zero():
    f=pd.read_csv(ROOT/'results/reciprocal_ANI.csv')
    assert (f.loc[f.status=='estimate_returned','ANI']>0).all()
    if (f.status=='no_estimate').any():assert f.loc[f.status=='no_estimate','ANI'].isna().all()
    for r in f.query("tool == 'skani 0.3.1'").itertuples():
        if r.status=='estimate_returned':assert 0<=r.Align_fraction_ref<=100 and 0<=r.Align_fraction_query<=100
def test_native_BOFdat_half_mass_is_exposed():
    r=json.loads((ROOT/'results/BOFdat_native_and_mass_correction.json').read_text())
    assert r['native_source_unmodified'];assert abs(r['native_fraction_of_requested_mass']-.5)<1e-10
    assert abs(r['corrected_mass_g_gDW']-r['DNA_fraction_g_gDW'])<1e-10
    c=r['correction_coefficients'];assert c['datp_c']==c['dttp_c'];assert c['dgtp_c']==c['dctp_c']
    assert abs(sum(c.values()))<1e-10
def test_candidate_changes_only_four_dNTPs_and_PPi():
    x=pd.read_csv(ROOT/'results/biomass_coefficients.csv').query("variant=='DNA_sequence_mass_corrected_candidate'")
    changed=set(x.loc[x.delta.abs()>1e-12,'metabolite'])
    assert changed=={'cpd00115_c0','cpd00357_c0','cpd00241_c0','cpd00356_c0','cpd00012_c0'}
def test_candidate_preserves_total_biomass_mass():
    f=pd.read_csv(ROOT/'results/biomass_candidate_phenotypes.csv').query("medium=='gut_like'").set_index('variant')
    assert abs(f.loc['baseline','net_biomass_mass_g_per_unit']-f.loc['DNA_sequence_mass_corrected_candidate','net_biomass_mass_g_per_unit'])<1e-10
    assert f.loc['DNA_native_raw_diagnostic','net_biomass_mass_g_per_unit']<.99
def test_precursor_producibility():
    f=pd.read_csv(ROOT/'results/biomass_precursor_tests.csv');assert len(f)>200;assert (f.solver_status=='optimal').all()
    blocked={name:set(x.loc[x.max_demand<1e-6,'metabolite']) for name,x in f.groupby('variant')}
    # Standalone drains for recycled carrier/polymer moieties are zero even in
    # the growth-capable baseline. The candidate must introduce no new failure.
    assert blocked['baseline']=={'cpd11493_c0','cpd00166_c0','cpd15665_c0'}
    assert blocked['DNA_sequence_mass_corrected_candidate']==blocked['baseline']
def test_no_positive_energy_dissipation_after_DNA_change():
    f=pd.read_csv(ROOT/'results/native_closed_boundary_energy_probes.csv');assert len(f)==24;assert (f.status=='optimal').all();assert f.maximum.abs().max()<1e-8
def test_native_GPR_knockouts_reproduce_v2():
    f=pd.read_csv(ROOT/'results/native_COBRA_GPR_knockouts.csv');assert f.difference.abs().max()<1e-7
    assert f.loc[f.design.isin(['OK2','OK3','OK4']),'native_growth'].abs().max()<1e-8
def test_native_reaction_knockouts_reproduce_v2():
    f=pd.read_csv(ROOT/'results/native_COBRA_reaction_knockouts.csv');assert f.difference.abs().max()<1e-7
    assert f.loc[f.design.isin(['OK2','OK3','OK4']),'native_growth'].min()>.6
def test_native_FVA_sensible_intervals():
    f=pd.read_csv(ROOT/'results/native_COBRA_exchange_FVA.csv');assert len(f)>100;assert (f.minimum<=f.maximum+1e-6).all()
def test_alignment_both_coverages_recorded():
    f=pd.read_csv(ROOT/'results/locus_evidence/alignment_evidence.csv');assert len(f)>500
    for c in ['query_coverage_percent','subject_coverage_percent']:assert f[c].between(0,100).all()
    assert json.loads((ROOT/'results/locus_evidence/status.json').read_text())['sequence_validated_mapping_to_numeric_model_genes']
def test_iME375_stereochemistry_not_conflated():
    f=pd.read_csv(ROOT/'results/comparators/medium_mapping_ledger.csv');x=f.query("model=='iME375' and ModelSEED=='cpd00159'")
    assert x.status.eq('missing_no_substitute').all()
    r=json.loads((ROOT/'results/comparators/structural_native_settings_audit.json').read_text());i=next(x for x in r if x['model']=='iME375')
    assert 'R047' in i['deposited_objective'];assert i['growth_objective']=='R617'
def test_historical_settings_are_recovered_with_source_cells():
    r=json.loads((ROOT/'results/original_StrainDesign_settings_recovered.json').read_text())
    assert r['recovered']['EPS_tilt']==.001;assert r['recovered']['candidate_count_recorded']==78;assert 13 in r['recovered']['source_cells']
def test_native_MICOM_LP_maxima_reproduce_v2():
    f=pd.read_csv(ROOT/'results/native_micom_linear/native_LP_maxima_vs_v2.csv')
    assert len(f)==60;assert f.status.eq('optimal').all();assert f.difference.abs().max()<1e-7
def test_native_MICOM_supplied_growth_states_and_accounting():
    f=pd.read_csv(ROOT/'results/native_micom_linear/native_fixed_growth_validation.csv')
    assert len(f)==48;assert f.status.eq('optimal').all()
    x=pd.read_csv(ROOT/'results/native_micom_linear/native_fixed_growth_exchange_accounting.csv')
    assert x.accounting_residual.abs().max()<1e-7
    for ctrl in ['block_both','no_lactate_block_both']:
        assert x.loc[(x.control==ctrl)&x.compound.isin(['cpd00159','cpd00221']),'B_specific'].max()<1e-7
