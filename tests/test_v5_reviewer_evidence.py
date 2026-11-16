from pathlib import Path
import json,hashlib,re
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]

def test_supplied_reports_all_text_preserved_and_every_block_answered():
    source=(ROOT/'inputs/reviewer_reports_user_supplied_2026-10-26.txt').read_text()
    record=json.loads((ROOT/'results/reviewer_v5/reviewer_comments_verbatim.json').read_text())
    replies=json.loads((ROOT/'results/reviewer_v5/point_by_point_responses.json').read_text())
    assert record['source_sha256']==hashlib.sha256(source.encode()).hexdigest()==replies['source_sha256']
    assert len(record['blocks'])==37
    assert {r['id'] for r in replies['blocks']}=={r['id'] for r in record['blocks']}
    mask=list(source)
    for block in record['blocks']:
        a,b=block['source_start_character'],block['source_end_character']
        assert source[a:b]==block['comment']
        assert hashlib.sha256(block['comment'].encode()).hexdigest()==block['comment_sha256']
        mask[a:b]=' '*len(block['comment'])
    leftovers={x.strip() for x in ''.join(mask).splitlines() if x.strip()}
    assert leftovers<={'Reviewer 1','Reviewer 2','Reviewer 3','General comment:','Specific comments:','Major comments','Minor comments','Comments','Attachment(s):'}

def test_new_reduced_subset_grows_and_all_targeted_OK2_OK4_genes_still_lethal():
    s=json.loads((ROOT/'results/reviewer_v5/inclusion_minimal_subset.json').read_text())
    assert s['additional_reactions_before']==237
    assert s['additional_reactions_retained']==112 and s['additional_reactions_disabled']==125
    assert s['native_growth']>.1 and s['mass_residual']<1e-8
    assert s['inclusion_minimal'] and not s['minimum_cardinality_claimed']
    cert=pd.read_csv(ROOT/'results/reviewer_v5/inclusion_minimal_subset_deletion_certificate.csv')
    assert len(cert)==112 and cert.growth_if_removed.max()<.1-1e-8
    d=pd.read_csv(ROOT/'results/reviewer_v5/inclusion_minimal_subset_global_GPR_retests.csv')
    assert d.loc[d.design=='WT','growth'].min()>.1
    assert d.loc[d.design.isin(['OK2','OK3','OK4']),'growth'].abs().max()<1e-7

def test_original_query_coverage_claim_has_pinned_source_and_actual_hits():
    d=pd.read_csv(ROOT/'results/reviewer_v5/original_gapseq_query_coverage_audit.csv')
    assert d.rows_with_coverage.tolist()==[20283,10613]
    assert d.min_recorded_query_coverage.eq(75).all() and d.rows_below_75.eq(0).all()
    for name in ['gapseq_find_c64db99b.sh','pinned_transporter.sh']:
        s=(ROOT/'inputs/reviewer_evidence_sources'/name).read_text()
        assert 'covcutoff=75' in s and '-qcov_hsp_perc $covcutoff' in s

def test_native_pathway_ranges_have_correct_native_multisubunit_GPR():
    d=pd.read_csv(ROOT/'results/reviewer_v5/native_Rnf_PFOR_acrylate_fluxes.csv')
    assert len(d)==21 and (d.minimum_90pct<=d.maximum_90pct+1e-7).all()
    assert d[d.reaction=='rxn05938_c0'].unique_gene_count.eq(0).all()
    x=d[d.reaction=='rxn02123_c0'];assert x.native_GPR.str.contains('and').all()

def test_partner_strain_verified_without_inventing_stereochemistry():
    d=json.loads((ROOT/'results/reviewer_v5/partner_literature_evidence.json').read_text())
    assert d['source_identical_to_deposited_partner'] and d['accession']=='GCF_000737885.1'
    assert not d['D_lactate_production_experiment_verified'] and not d['net_L_lactate_consumption_experiment_verified']

def test_current_index_has_one_definition_for_every_S1_to_S73_and_real_targets():
    d=pd.read_csv(ROOT/'tables/SUPPLEMENTARY_INDEX_v5.csv')
    assert len(d)==73 and set(d.source_id)=={'S'+str(i) for i in range(1,74)}
    for row in d.itertuples():
        for p in row.paths.split(';'):assert (ROOT/p).exists(),(row.source_id,p)

def test_complete_report_recovery_is_resolved_but_unperformed_requests_are_not():
    d=json.loads((ROOT/'results/reviewer_v5/status_summary.json').read_text())
    assert d['missing_report_text_resolved'] and not d['all_requests_fully_satisfied']
    r=json.loads((ROOT/'results/reviewer_v5/point_by_point_responses.json').read_text())['blocks']
    assert next(x for x in r if x['id']=='R3.M2')['status']=='unresolved_full_NM10_BL7_simulations'

def test_response_line_references_match_actual_clean_aux_and_all_comments_render():
    r=json.loads((ROOT/'results/reviewer_v5/reviewer_coverage_with_actual_lines.json').read_text())['blocks']
    aux=(ROOT/'manuscript/manuscript_clean_v5/main.aux').read_text()
    labels={k:int(v) for k,v in re.findall(r'\\newlabel\{(v[2345]-[^}]+)\}\{\{\\getpagewiselinenumber\s*\{(\d+)\}',aux)}
    assert len(r)==37
    for block in r:
        assert len(block['compiled_clean_line_ranges'])==len(block['line_label_prefixes'])
        for ref in block['compiled_clean_line_ranges']:
            assert ref['start']==labels[ref['label_prefix']+'-start'] and ref['end']==labels[ref['label_prefix']+'-end']
    s=(ROOT/'manuscript/response_v5/response_letter_npjSB_v5.tex').read_text()
    assert s.count(r'\RevPoint{')==37 and r'\VFiveLines' not in s

def test_final_table_numbers_are_consecutive_and_not_double_counted_by_ltablex():
    p=ROOT/'manuscript/manuscript_clean_v5'
    assert 'ltablex' not in (p/'main.tex').read_text()
    labels={k:int(v) for k,v in re.findall(r'\\newlabel\{(tab:[^}]+)\}\{\{(\d+)\}',(p/'main.aux').read_text())}
    assert len(labels)==7 and sorted(labels.values())==list(range(1,8))
    assert labels['tab:strain_design']==6
