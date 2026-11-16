from common import *
import pandas as pd,nbformat,re

def main():
    index=pd.read_csv(ROOT/'tables/SUPPLEMENTARY_INDEX_v6.csv');assert len(index)==87 and index.source_id.nunique()==87
    for paths in index.paths:
        for p in paths.split(';'):assert (ROOT/p).is_file(),p
    rep=json.loads((ROOT/'results/reviewer_v6/reviewer_coverage_with_actual_lines.json').read_text());source=(ROOT/'inputs/reviewer_reports_user_supplied_2026-10-26.txt').read_text()
    assert sha(ROOT/'inputs/reviewer_reports_user_supplied_2026-10-26.txt')==rep['source_sha256'];assert len(rep['blocks'])==37
    for b in rep['blocks']:
        assert source[b['source_start_character']:b['source_end_character']]==b['comment']
        assert all(x['end']>=x['start']>0 for x in b['compiled_clean_line_ranges'])
    nb=nbformat.read(ROOT/'notebooks/Review_v6_EXECUTED.ipynb',as_version=4)
    assert sum(c.cell_type=='code' for c in nb.cells)==7
    assert all(o.output_type!='error' for c in nb.cells if c.cell_type=='code' for o in c.outputs)
    # All inherited scientific inputs/results are byte-identical to frozen v5.
    old=pd.read_csv(ROOT/'provenance/inherited_MANIFEST_v5.tsv',sep='\t');checked=0;failures=[]
    for row in old.itertuples():
        if not row.path.startswith(('baseline_v2/','inputs/','models/','results/','figures/')):continue
        p=ROOT/row.path
        if not p.exists() or sha(p)!=row.sha256:failures.append(row.path)
        checked+=1
    assert not failures,failures
    figs=[]
    for p in (ROOT/'manuscript/manuscript_clean_v5').iterdir():
        if re.match(r'Figure_[3-7]',p.name) and p.is_file():
            q=ROOT/'manuscript/manuscript_clean_v6'/p.name;assert sha(p)==sha(q);figs.append(p.name)
    settings=[]
    for p in ROOT.glob('results/matched_micom_v6/*/settings.json'):
        x=json.loads(p.read_text());assert sha(ROOT/x['M'])==x['M_sha256'] and sha(ROOT/x['B'])==x['B_sha256'];settings.append(str(p.relative_to(ROOT)))
    assert len(settings)==8
    record('provenance/v6_release_integrity_checks.json',{'verified_utc':now(),'supplementary_sources':87,'reviewer_blocks':37,'notebook_code_cells':7,'notebook_errors':0,'inherited_scientific_files_byte_verified':checked,'inherited_scientific_mismatches':failures,'Figures_3_7_preserved':figs,'native_community_input_hashes_verified':settings,'passed':True})
    print('Release source integrity passed;',checked,'inherited scientific files unchanged;87 source entries;37 comments;8 native input sets.')
if __name__=='__main__':main()
