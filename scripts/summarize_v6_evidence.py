from common import *
import pandas as pd,numpy as np,shutil

def aggregate(pattern,dest):
    fs=sorted(ROOT.glob(pattern));assert fs,pattern
    data=pd.concat([pd.read_csv(f).assign(source_file=str(f.relative_to(ROOT))) for f in fs],ignore_index=True)
    data.to_csv(ROOT/dest,index=False);return data

def main():
    out=ROOT/'results/reviewer_v6';out.mkdir(exist_ok=True)
    mono=aggregate('results/matched_monocultures_v6/*/phenotypes.csv','tables/S75_matched_monoculture_phenotypes_v6.csv')
    aa=aggregate('results/matched_monocultures_v6/*/AA_singles_and_pairs.csv','tables/S76_matched_AA_singles_pairs_v6.csv')
    grid=aggregate('results/matched_monocultures_v6/*/maintenance_uptake_grid.csv','tables/S77_matched_maintenance_grid_v6.csv')
    designs=aggregate('results/matched_monocultures_v6/*/intervention_implementations.csv','tables/S78_matched_interventions_v6.csv')
    fva=aggregate('results/matched_monocultures_v6/*/exchange_FVA.csv','tables/S79_matched_exchange_FVA_v6.csv')
    native=aggregate('results/matched_micom_v6/*/native_tradeoff_controls.csv','tables/S80_matched_native_MICOM_v6.csv')
    cfva=aggregate('results/matched_micom_v6/*/exchange_FVA95.csv','tables/S81_matched_community_FVA_v6.csv')
    shutil.copy2(ROOT/'results/functional_curation_v6/GPR_edit_ledger.csv',ROOT/'tables/S82_primary_GPR_curation_v6.csv')
    shutil.copy2(ROOT/'results/functional_curation_v6/B22L_alignment_evidence.csv',ROOT/'tables/S83_22L_sequence_evidence_v6.csv')
    shutil.copy2(ROOT/'results/matched_biomass_v6/growth_maintenance_sensitivity.csv',ROOT/'tables/S84_matched_biomass_sensitivity_v6.csv')
    shutil.copy2(ROOT/'results/comparator_calibration_v6/common_carbon_budget_simulations.csv',ROOT/'tables/S85_comparator_calibration_v6.csv')
    shutil.copy2(ROOT/'results/comparator_calibration_v6/comparator_ATP_cycle_checks.csv',ROOT/'tables/S86_comparator_energy_checks_v6.csv')
    audits=[json.loads(p.read_text()) for p in sorted((ROOT/'results/matched_monocultures_v6').glob('*/model_audit.json'))]
    pd.DataFrame(audits).to_csv(ROOT/'tables/S74_new_matched_reconstruction_audit_v6.csv',index=False)
    shown=mono[mono.scenario.isin(['gut_like','carbon_rich','AA_rich'])][['model','scenario','growth','cpd00211','cpd00141','cpd00159','cpd00221']]
    shown.to_csv(out/'matched_media_summary.csv',index=False)
    selected=native[native.M_abundance.eq(.5)&native.tradeoff.eq(.95)&native.control.isin(['baseline','no_external_lactate','no_lactate_block_both'])]
    selected.to_csv(out/'community_summary.csv',index=False)
    selected_flux=[];certs=[];account=[]
    for d in sorted((ROOT/'results/matched_micom_v6').iterdir()):
        if not d.is_dir():continue
        f=pd.read_csv(d/'native_taxon_community_fluxes.csv');f=f[f.M_abundance.eq(.5)&f.tradeoff.eq(.95)&f.control.isin(['baseline','no_external_lactate','no_lactate_block_both'])]
        selected_flux.append(f)
        c=pd.read_csv(d/'QP_certificates.csv');c['pair']=d.name;certs.append(c)
        allf=pd.read_csv(d/'native_taxon_community_fluxes.csv');account.append({'pair':d.name,'max_pool_residual':float(allf.accounting_residual.abs().max())})
    pd.concat(selected_flux).to_csv(out/'selected_community_fluxes.csv',index=False)
    cert=pd.concat(certs);cert.to_csv(out/'all_QP_certificates.csv',index=False)
    before=native[native.pair.eq('GPRcurated_original22L')];after=native[native.pair.eq('GPRcurated_filtered22L')]
    comp=before.merge(after,on=['M_abundance','tradeoff','control'],suffixes=('_before','_after'))
    maxdiff={c:float((comp[c+'_after']-comp[c+'_before']).abs().max()) for c in ['community_max','community_growth','growth_M','growth_B']}
    facts={'generated_utc':now(),'new_reconstructions':audits,'mono_rows':len(mono),'AA_rows':len(aa),'maintenance_grid_rows':len(grid),'design_rows':len(designs),
           'ordinary_exchange_FVA_rows':len(fva),'native_community_cases':len(native),'native_community_optimal_cases':int(native.status.eq('optimal').sum()),'community_FVA_rows':len(cfva),
           'native_max_mass_residual':float(native.mass_residual.max()),'native_max_bound_violation':float(native.bound_violation.max()),
           'QP_max_global_gap':float(native.QP_global_gap.max()),'QP_max_stationarity':float(native.QP_stationarity.max()),
           'pool_accounting':account,'22L_DLDH_filter_max_growth_changes':maxdiff,
           'experimental_MJR_biomass_available':False,'experimental_22L_enantiomer_assay_verified':False,
           'new_models_are_published_NM10_BL7_models':False,'NM10_Bcd':'Public misc_feature at APHY01000022.1:15273-16423 minus strand explicitly notes potential frameshift; new ORFs are partial high-identity Bcd homologs. Strict Bcd-disabled sensitivity, not proof of biological loss.'}
    record('results/reviewer_v6/facts.json',facts)
    statuses=[{'task':'NM10 and BL7 matched comparison','status':'EXECUTED_NEW_RECONSTRUCTIONS','actual_result':'Uniform gapseq2.1.0/DB1.5 models, all media and full matched carbon panel, 20-AA singles/190 pairs, maintenance and intervention controls, native communities; published old ModelSEED models remain unavailable.'},
              {'task':'Functional pathway curation','status':'EXECUTED_SUPPORTED_GPR_EDITS_AND_UNCERTAINTY_CONTROLS','actual_result':'Bcd-Etf and six-subunit Rnf rules corrected; PFOR GPRs inferred from full-length references. Acrylate catalytic specificity, complete Lar maturation and terminal C5/C6 chemistry remain unresolved; no invented additions.'},
              {'task':'22L stereochemistry','status':'GENOME_AUDITED_AND_MODEL_ARTIFACT_TESTED','actual_result':'All349 deposited gene accession versions verified; strong L-LDH evidence; unsupported D-LDH gapfill disabled in sensitivity; native community effects calculated. No strain-matched enantiomer-resolved experiment found.'},
              {'task':'Biomass and comparator calibration','status':'EXECUTED_BOUNDED_CALIBRATION','actual_result':'Sequence-informed biomass candidates, mass/precursor/ATP tests; four AGORA formula-mass candidates; common rich supplies and gross organic-C budget; iME375 source inconsistencies and ATP cycle corrected in separate candidates. Missing target measurements and generic-polymer masses prevent empirical full calibration.'},
              {'task':'Submission readiness','status':'COMPUTATIONAL_REVISION_FOR_AUTHOR_REVIEW','actual_result':'New numerical results replace formerly unperformed comparisons. Experimental validation and residual pathway specificity are limitations; no unconditional submission-ready claim.'}]
    pd.DataFrame(statuses).to_csv(ROOT/'TASK_STATUS_v6.csv',index=False)
    print(json.dumps({k:v for k,v in facts.items() if k not in ['new_reconstructions','pool_accounting']},indent=2))

if __name__=='__main__':main()
