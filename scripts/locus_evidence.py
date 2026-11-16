"""Curated-reference homology candidates; never equate a hit with activity."""
from common import *
import json,re,collections
import pandas as pd
from Bio import SeqIO

def main():
    refs=ROOT/'inputs/reference_proteins/combined_reference_panel.faa';proteins=ROOT/'inputs/genomes/GCA_001546855.1/protein.faa'
    out=ROOT/'results/locus_evidence';out.mkdir(exist_ok=True)
    run(['diamond','makedb','--in',proteins,'--db',ROOT/'inputs/target_proteins'],'alignments/target_makedb')
    run(['diamond','makedb','--in',refs,'--db',ROOT/'inputs/reviewed_refs'],'alignments/references_makedb')
    fields=['qseqid','sseqid','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore','qseq','sseq']
    hit=out/'reference_to_MJR8396C.tsv';reverse=out/'MJR8396C_to_reference.tsv'
    run(['diamond','blastp','--query',refs,'--db',ROOT/'inputs/target_proteins','--sensitive','--threads','4','--evalue','1e-5','--max-target-seqs','5','--outfmt','6',*fields,'--out',hit],'alignments/forward',300)
    run(['diamond','blastp','--query',proteins,'--db',ROOT/'inputs/reviewed_refs','--sensitive','--threads','4','--evalue','1e-5','--max-target-seqs','1','--outfmt','6',*fields,'--out',reverse],'alignments/reverse',300)
    h=pd.read_csv(hit,sep='\t',names=fields);back=pd.read_csv(reverse,sep='\t',names=fields)
    def acc(s):return s.split('|')[1] if '|' in s else s
    h['reference_accession']=h.qseqid.map(acc);back['reference_accession']=back.sseqid.map(acc)
    rb=dict(zip(back.qseqid,back.reference_accession));h['reciprocal_best_hit']=[rb.get(s)==q for s,q in zip(h.sseqid,h.reference_accession)]
    h['query_coverage_percent']=100*(h.qend-h.qstart+1)/h.qlen;h['subject_coverage_percent']=100*(h.send-h.sstart+1)/h.slen
    data=json.loads((ROOT/'inputs/reference_proteins/combined_reference_panel.json').read_text())['results'];metadata={x['primaryAccession']:x for x in data};refrows=[]
    for a,x in metadata.items():
        desc=x['proteinDescription'];name=desc.get('recommendedName',{}).get('fullName',{}).get('value','')
        genes=';'.join(g.get('geneName',{}).get('value','') for g in x.get('genes',[]))
        ecs=';'.join(e['value'] for e in desc.get('recommendedName',{}).get('ecNumbers',[]))
        pubs=[db['id'] for r in x.get('references',[]) for db in r.get('citation',{}).get('citationCrossReferences',[]) if db['database']=='PubMed']
        domains=';'.join(d['id'] for d in x.get('uniProtKBCrossReferences',[]) if d['database'] in ['Pfam','InterPro'])
        refrows.append({'reference_accession':a,'sequence_version':x.get('entryAudit',{}).get('sequenceVersion'),'entry_version':x.get('entryAudit',{}).get('entryVersion'),'reviewed':x.get('entryType'),'protein_name':name,'reference_genes':genes,'EC':ecs,'reference_organism':x.get('organism',{}).get('scientificName'),'PMIDs':';'.join(pubs),'reference_domains':domains,'protein_existence':x.get('proteinExistence'),'source_url':'https://www.uniprot.org/uniprotkb/'+a})
    meta=pd.DataFrame(refrows);meta.to_csv(out/'reference_provenance.csv',index=False);h=h.merge(meta,on='reference_accession',how='left')
    loc=pd.read_csv(ROOT/'results/genbank_CDS_locus_ledger.csv').query("assembly == 'GCA_001546855.1'")
    h=h.merge(loc,left_on='sseqid',right_on='protein_id',how='left')
    h['candidate_classification']=h.apply(lambda r:'homology_candidate' if r.pident>=30 and r.query_coverage_percent>=70 and r.subject_coverage_percent>=70 else 'partial_or_distant_hit_requires_review',axis=1)
    h['curation_decision']='No reaction/direction/GPR edit: homology is not activity; mapping to original numeric JGI genes remains unresolved.'
    catalytic=[]
    for r in h.itertuples():
        x=metadata[r.reference_accession];features=[f for f in x.get('features',[]) if f['type'] in ['Active site','Binding site']]
        mapping={};qi=int(r.qstart);si=int(r.sstart)
        for q,s in zip(r.qseq,r.sseq):
            if q!='-':mapping[qi]=(s,si if s!='-' else None,q);qi+=1
            if s!='-':si+=1
        for feat in features:
            lo=feat['location']['start'].get('value');hi=feat['location']['end'].get('value')
            if lo!=hi or lo not in mapping:continue
            residue,pos,original=mapping[lo]
            catalytic.append({'reference_accession':r.reference_accession,'protein_id':r.sseqid,'type':feat['type'],'reference_position':lo,'reference_residue':original,'target_position':pos,'target_residue':residue,'conserved':original==residue,'description':feat.get('description',''),'evidence':json.dumps(feat.get('evidences',[]))})
    h.drop(columns=['qseq','sseq']).to_csv(out/'alignment_evidence.csv',index=False)
    pd.DataFrame(catalytic).to_csv(out/'reference_annotated_residue_comparisons.csv',index=False)
    good=h.query("candidate_classification == 'homology_candidate'").copy()
    # Explicit local gene neighborhood for every well-covered candidate.
    near=[]
    for r in good.drop_duplicates('protein_id').itertuples():
        for n in loc[(loc.scaffold==r.scaffold)&(loc.start_1based<=r.end_1based_inclusive+10000)&(loc.end_1based_inclusive>=r.start_1based-10000)].itertuples():
            near.append({'candidate_protein':r.protein_id,'candidate_locus':r.locus_tag,'neighbor_locus':n.locus_tag,'neighbor_protein':n.protein_id,'scaffold':n.scaffold,'start':n.start_1based,'end':n.end_1based_inclusive,'strand':n.strand,'product':n.product})
    pd.DataFrame(near).to_csv(out/'candidate_neighborhoods_10kb.csv',index=False)
    # Preserve all pre-existing model annotations, explicitly as proxy evidence.
    import cobra
    cobra.Configuration().solver='glpk';m=cobra.io.read_sbml_model(str(BASE/'models/Megasphaera_sp_MJR8396C_filled_anaerobic.xml'))
    rows=[]
    for r in m.reactions:
        if not re.search('lactate|lactoyl|acrylo|butyryl|ferredoxin|pyruvate.*oxid|rnf|electron.transfer|croton|CoA.transfer',r.name,re.I):continue
        rows.append({'reaction':r.id,'name':r.name,'current_GPR':r.gene_reaction_rule,'bounds':str(r.bounds),'reaction_stoichiometry':r.reaction,'model_gene_annotations':json.dumps({g.id:g.annotation for g in r.genes}),'public_locus_mapping':'unresolved_without_original_JGI_protein_sequences','alignment_table':'alignment_evidence.csv','curation_decision':'preserve original stoichiometry/direction/GPR; candidates in separate evidence table, no verified activity claimed'})
    pd.DataFrame(rows).to_csv(out/'reaction_GPR_curation_ledger.csv',index=False)
    # Completeness is only a candidate inventory, never a complex assignment.
    inventory=[]
    for name,regex in [('L_LDH','L-lactate dehydrogenase'),('D_LDH','D-lactate dehydrogenase'),('LarA','lactate racemase'),('LarB','larB'),('LarC','larC'),('LarE','larE'),('Bcd','butyryl-CoA dehydrogenase'),('EtfA','etfA'),('EtfB','etfB'),('PFOR','pyruvate.*(?:ferredoxin|oxidoreductase)'),('RnfA','rnfA'),('RnfB','rnfB'),('RnfC','rnfC'),('RnfD','rnfD'),('RnfE','rnfE'),('RnfG','rnfG'),('lactoyl_CoA_dehydratase','lactoyl-CoA dehydratase'),('CoA_transferase','CoA.transferase')]:
        subset=good[(good.protein_name.fillna('')+' '+good.reference_genes.fillna('')).str.contains(regex,case=False,regex=True)]
        inventory.append({'component':name,'candidate_loci':';'.join(sorted(set(subset.locus_tag.dropna()))),'candidate_count':subset.locus_tag.nunique(),'reference_hits':';'.join(sorted(set(subset.reference_accession))),'interpretation':'Candidates only; absence from this curated reference panel does not prove genomic absence; no full complex assigned.'})
    pd.DataFrame(inventory).to_csv(out/'complex_component_candidate_inventory.csv',index=False)
    record('results/locus_evidence/status.json',{'reference_sequences':len(data),'alignments':len(h),'well_covered_candidates':good.protein_id.nunique(),'sequence_validated_mapping_to_numeric_model_genes':False,'activity_verified':False,'models_edited':False,'limitations':['Curated UniProt panel is not an exhaustive domain profile search.','Reciprocal best hit is relative to the selected reference panel.','Required Lar maturation, substrate specificity and multifunctional-complex completeness require expert review.','No original JGI proteins were found in approved project folders.']})
if __name__=='__main__':main()
