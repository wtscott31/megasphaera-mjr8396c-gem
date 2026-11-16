"""Apply the same sequence and complex logic to all NEW matched models."""
from functional_curation_v6 import *
from matched_monocultures_v6 import minimal_off,prune
import gzip,itertools

OUT=ROOT/'results/matched_functional_curation_v6'

def main():
    OUT.mkdir(exist_ok=True);ref=ROOT/'inputs/reference_proteins/combined_reference_panel.faa'
    meta=pd.read_csv(ROOT/'results/locus_evidence/reference_provenance.csv');all_edits=[];all_hits=[];deletions=[]
    fields=['qseqid','sseqid','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore']
    for n in ['MJR8396C','NM10','BL7','B22L']:
        d=ROOT/'results/matched_reconstruction_v6'/n
        fasta=OUT/(n+'.faa')
        with gzip.open(d/(n+'.faa.gz'),'rt') as f:fasta.write_text(f.read())
        records=list(SeqIO.parse(fasta,'fasta'));coords={}
        for rec in records:
            desc=rec.description.split('#')
            coords[rec.id]={'contig':rec.id.rsplit('_',1)[0],'start':int(desc[1]),'end':int(desc[2]),'strand':int(desc[3])}
        run(['diamond','makedb','--in',fasta,'--db',OUT/n],'v6_matched_'+n+'_db')
        run(['diamond','blastp','--query',ref,'--db',OUT/n,'--sensitive','--threads','4','--evalue','1e-5','--max-target-seqs','5','--outfmt','6',*fields,'--out',OUT/(n+'_reference_hits.tsv')],'v6_matched_'+n+'_references')
        h=pd.read_csv(OUT/(n+'_reference_hits.tsv'),sep='\t',names=fields);h['reference_accession']=h.qseqid.str.split('|').str[1]
        h['query_coverage_percent']=100*(h.qend-h.qstart+1)/h.qlen;h['subject_coverage_percent']=100*(h.send-h.sstart+1)/h.slen
        h=h.merge(meta,on='reference_accession',how='left');h['model']=n
        h['model_gene']=h.sseqid.str.replace('.','_',regex=False);all_hits.append(h)
        m=load(path=d/(n+'.xml'));before=m.copy()
        good=h[(h.pident>=30)&(h.query_coverage_percent>=75)&(h.subject_coverage_percent>=75)].sort_values('bitscore',ascending=False)
        edits=[]
        if n=='B22L':
            hits=good[good.EC.fillna('').str.split(';').map(lambda x:'1.1.1.28' in x)]
            assert len(hits)==0 and not m.reactions.rxn00500_c0.gene_reaction_rule
            edits.append({'reaction':'rxn00500_c0','before_GPR':'','after_GPR':'','before_bounds':str(m.reactions.rxn00500_c0.bounds),'after_bounds':'(0, 0)','reason':'Unsupported D-NAD-LDH gapfill disabled in a separate sensitivity; no well-covered curated reference candidate. Not proof of genomic absence.'})
            m.reactions.rxn00500_c0.bounds=(0,0)
        else:
            # Select coherent Bcd/Etf neighborhood using a characterized Bcd.
            complete=good[good.reference_accession.eq('Q06319')&(good.pident>=50)&(good.query_coverage_percent>=90)&(good.subject_coverage_percent>=90)]
            complex_genes=[];complex_loci=[]
            for bcd in complete.itertuples():
                bc=coords[bcd.sseqid];cg=[bcd.model_gene];cl=[bcd.sseqid]
                for component in ['etfA','etfB']:
                    candidates=good[good.reference_genes.fillna('').eq(component)].drop_duplicates('sseqid')
                    candidates=candidates[candidates.sseqid.map(lambda s:coords[s]['contig']==bc['contig'] and abs(coords[s]['start']-bc['start'])<5000)]
                    if len(candidates)==0:break
                    chosen=candidates.iloc[0];cg.append(chosen.model_gene);cl.append(chosen.sseqid)
                if len(cg)==3:complex_genes=cg;complex_loci=cl;break
            r=m.reactions.rxn90001_c0;old=r.gene_reaction_rule
            if complex_genes:
                r.gene_reaction_rule=' and '.join(complex_genes)
                edits.append({'reaction':r.id,'before_GPR':old,'after_GPR':r.gene_reaction_rule,'support_proteins':';'.join(complex_loci),'reason':'Complete Bcd-Etf alpha/beta neighborhood; catalytic Bcd required, ETF-only OR rule corrected.'})
            else:
                # A split ORF is not presumed to make an active enzyme. Remote
                # acyl-CoA-dehydrogenase homologs may act in trans; unresolved.
                partial=h[h.reference_accession.eq('Q06319')&(h.pident>=90)]
                partial.to_csv(OUT/(n+'_Bcd_high_identity_fragments.csv'),index=False)
                r.gene_reaction_rule='';r.bounds=(0,0)
                edits.append({'reaction':r.id,'before_GPR':old,'after_GPR':'','before_bounds':str(before.reactions.rxn90001_c0.bounds),'after_bounds':'(0, 0)',
                              'support_proteins':';'.join(sorted(set(partial.sseqid))),
                              'reason':'STRICT SENSITIVITY ONLY: no full-length Bcd homolog in the Etf neighborhood; high-identity split ORFs recorded. Bcd blocked rather than assuming fragment complementation. Remote homolog substrate activity and assembly error remain possible; genomic absence NOT claimed.'})
            # A complete Rnf operon must contain all six homologous components.
            cn=good[good.reference_genes.fillna('').eq('rnfC')].iloc[0];cc=coords[cn.sseqid];rg=[];rp=[]
            for component in ['rnfC','rnfD','rnfG','rnfE','rnfA','rnfB']:
                candidates=good[good.reference_genes.fillna('').eq(component)].drop_duplicates('sseqid')
                candidates=candidates[candidates.sseqid.map(lambda s:coords[s]['contig']==cc['contig'] and abs(coords[s]['start']-cc['start'])<10000)]
                assert len(candidates)>0,(n,component)
                chosen=candidates.iloc[0];rg.append(chosen.model_gene);rp.append(chosen.sseqid)
            r=m.reactions.rxn90053_c0;old=r.gene_reaction_rule;r.gene_reaction_rule=' and '.join(rg)
            edits.append({'reaction':r.id,'before_GPR':old,'after_GPR':r.gene_reaction_rule,'support_proteins':';'.join(rp),'reason':'Six-component RnfCDGEAB AND rule; same conserved-neighborhood test for each strain.'})
            pfor=good[good.reference_accession.isin(['Q2RMD6','P94692'])&(good.pident>=50)&(good.query_coverage_percent>=90)&(good.subject_coverage_percent>=90)]
            genes=sorted(set(pfor.model_gene));assert genes
            r=m.reactions.rxn05938_c0;old=r.gene_reaction_rule;r.gene_reaction_rule=' or '.join(genes)
            edits.append({'reaction':r.id,'before_GPR':old,'after_GPR':r.gene_reaction_rule,'support_proteins':';'.join(sorted(set(pfor.sseqid))),'reason':'Full-length PFOR homologs; inferred isozymes, not experimentally verified substrate activities.'})
        dest=ROOT/'models'/(n+'_matched_sequence_curated_v6.xml');m.id=n+'_matched_sequence_curated_v6'
        m.notes['scope']='NEW matched gapseq reconstruction, followed by the same explicit curated-reference/complex rules. Gene inference is not an activity assay. B22L D-LDH restriction is a conservative evidence-filtered sensitivity.'
        cobra.io.write_sbml_model(m,str(dest))
        for e in edits:all_edits.append({'model':n,'source_sha256':sha(d/(n+'.xml')),'candidate_sha256':sha(dest),**e})
        for source,mod in [('uncurated',before),('curated',m)]:
            medium(mod)
            for rid in ['rxn90001_c0','rxn90053_c0','rxn05938_c0','rxn00916_c0','rxn00410_c0']:
                if rid not in mod.reactions:continue
                r=mod.reactions.get_by_id(rid)
                if not r.gene_reaction_rule:
                    deletions.append({'model':n,'source':source,'reaction':rid,'status':'unassigned_GPR'});continue
                for gs in minimal_off(r.gpr.body):
                    with mod:
                        cobra.manipulation.knock_out_model_genes(mod,list(gs))
                        deletions.append({'model':n,'source':source,'reaction':rid,'genes':';'.join(sorted(gs)),'reaction_disabled':mod.reactions.get_by_id(rid).bounds==(0,0),**phenotype(mod)})
        print(n,'matched sequence curation complete',flush=True)
    pd.DataFrame(all_edits).to_csv(OUT/'matched_GPR_and_bound_edits.csv',index=False)
    pd.concat(all_hits).to_csv(OUT/'matched_reference_alignment_evidence.csv',index=False)
    pd.DataFrame(deletions).to_csv(OUT/'matched_GPR_essentiality_checks.csv',index=False)

if __name__=='__main__':main()
