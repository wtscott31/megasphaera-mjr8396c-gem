"""Sequence-supported GPR corrections and explicit 22L evidence filtering.

Homology supports model annotation, not experimental enzyme validation.
No new racemase, chain-elongation chemistry, or kinetic bounds are invented.
"""
from native_analysis import *
from locus_evidence import main as unused_locus_main

OUT=ROOT/'results/functional_curation_v6'
M_PATH=ROOT/'models/MJR8396C_sequence_supported_GPR_curated_v6.xml'
B_PATH=ROOT/'models/B22L_D_LDH_unsupported_gapfill_disabled_v6.xml'

def b22l_evidence():
    refs=ROOT/'inputs/reference_proteins/combined_reference_panel.faa'
    proteins=ROOT/'inputs/genomes/GCF_000737885.1/protein.faa'
    run(['diamond','makedb','--in',proteins,'--db',OUT/'B22L'],'v6_22L_makedb')
    fields=['qseqid','sseqid','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore','qseq','sseq']
    run(['diamond','blastp','--query',refs,'--db',OUT/'B22L','--sensitive','--threads','4','--evalue','1e-5','--max-target-seqs','5','--outfmt','6',*fields,'--out',OUT/'reference_to_B22L.tsv'],'v6_22L_homology')
    run(['diamond','blastp','--query',proteins,'--db',ROOT/'inputs/reviewed_refs','--sensitive','--threads','4','--evalue','1e-5','--max-target-seqs','1','--outfmt','6',*fields,'--out',OUT/'B22L_to_reference.tsv'],'v6_22L_reverse')
    h=pd.read_csv(OUT/'reference_to_B22L.tsv',sep='\t',names=fields)
    back=pd.read_csv(OUT/'B22L_to_reference.tsv',sep='\t',names=fields)
    def acc(x):return x.split('|')[1] if '|' in x else x
    h['reference_accession']=h.qseqid.map(acc);rb=dict(zip(back.qseqid,back.sseqid.map(acc)))
    h['reciprocal_best_within_panel']=[rb.get(s)==q for s,q in zip(h.sseqid,h.reference_accession)]
    h['query_coverage_percent']=100*(h.qend-h.qstart+1)/h.qlen
    h['subject_coverage_percent']=100*(h.send-h.sstart+1)/h.slen
    meta=pd.read_csv(ROOT/'results/locus_evidence/reference_provenance.csv')
    h=h.merge(meta,on='reference_accession',how='left');loc=[]
    for rec in SeqIO.parse(ROOT/'inputs/genomes/GCF_000737885.1/genomic.gbff','genbank'):
        for f in rec.features:
            q=f.qualifiers
            if f.type!='CDS' or 'protein_id' not in q:continue
            loc.append({'protein_id':q['protein_id'][0],'locus_tag':q.get('locus_tag',[''])[0],
                        'product':q.get('product',[''])[0],'scaffold':rec.id,
                        'start_1based':int(f.location.start)+1,'end_1based':int(f.location.end),
                        'strand':f.location.strand,'translation_length':len(q.get('translation',[''])[0])})
    loc=pd.DataFrame(loc);loc.to_csv(OUT/'B22L_CDS_ledger.csv',index=False)
    h=h.merge(loc,left_on='sseqid',right_on='protein_id',how='left')
    h['well_covered_candidate']=(h.pident>=30)&(h.query_coverage_percent>=75)&(h.subject_coverage_percent>=75)
    h.drop(columns=['qseq','sseq']).to_csv(OUT/'B22L_alignment_evidence.csv',index=False)
    b=load('B');ps={p.id:str(p.seq) for p in SeqIO.parse(proteins,'fasta')}
    rows=[]
    for g in b.genes:
        versioned=re.sub(r'_(\d+)$',r'.\1',g.id)
        rows.append({'model_gene':g.id,'NCBI_protein':versioned,'present_exact_accession_version':versioned in ps,
                     'AA_sha256':hashlib.sha256(ps[versioned].encode()).hexdigest() if versioned in ps else None,
                     'locus_tags':';'.join(loc.loc[loc.protein_id.eq(versioned),'locus_tag'])})
    pd.DataFrame(rows).to_csv(OUT/'B22L_deposited_gene_accession_crosscheck.csv',index=False)
    inventory=[]
    for name,mask in [('L_NAD_LDH',h.EC.fillna('').str.contains(r'(^|;)1\.1\.1\.27($|;)',regex=True)),
                      ('D_NAD_LDH',h.EC.fillna('').str.contains(r'(^|;)1\.1\.1\.28($|;)',regex=True)),
                      ('LarA',h.reference_genes.fillna('').eq('larA')),
                      ('LarB',h.reference_genes.fillna('').eq('larB')),
                      ('LarC',h.reference_genes.fillna('').eq('larC')),
                      ('LarE',h.reference_genes.fillna('').eq('larE'))]:
        subset=h[mask&h.well_covered_candidate].sort_values('bitscore',ascending=False)
        inventory.append({'component':name,'candidate_loci':';'.join(sorted(set(subset.locus_tag.dropna()))),
                          'candidate_proteins':';'.join(sorted(set(subset.sseqid))),
                          'best_hit':subset.drop(columns=['qseq','sseq']).head(1).to_dict('records'),
                          'meaning':'Limited curated-panel homology evidence; non-detection does not prove genomic absence'})
    record('results/functional_curation_v6/B22L_inventory.json',inventory)
    # Isolate an unsupported gapfill rather than assigning the L-LDH to D-LDH.
    original=b.reactions.rxn00500_c0
    assert not original.gene_reaction_rule
    bounds=original.bounds;original.bounds=(0,0)
    b.id='B22L_evidence_filtered_D_LDH_gapfill_disabled_v6'
    b.notes['scope']='Sensitivity model. D-NAD-LDH rxn00500 disabled: deposited reaction has no GPR; no well-covered D-NAD-LDH candidate in curated panel. This is not experimental absence evidence. L-NAD-LDH remains reversible; no unsupported physiological direction imposed.'
    cobra.io.write_sbml_model(b,str(B_PATH))
    record('results/functional_curation_v6/B22L_curation_decision.json',{'reaction':'rxn00500_c0','original_bounds':bounds,'candidate_bounds':[0,0],
           'source_model_sha256':sha(FILES['B']),'current_annotation_proteins':len(ps),'deposited_genes':len(rows),
           'present_accession_versions':sum(x['present_exact_accession_version'] for x in rows),
           'experimental_stereochemistry_verified':False,
           'decision':'Keep deposited model; release evidence-filtered sensitivity. D secretion and L uptake are not established 22L traits.',
           '22L_primary_fermentation_source':'https://doi.org/10.3389/fmicb.2015.01030',
           'assay_limitation':'Total lactic-acid HPLC with RI detector, not enantiomer-resolved. No verified 22L net L-lactate uptake assay located.'})

def m_curation():
    m=load(path=ROOT/'models/MJR8396C_sequence_annotated_primary.xml')
    original=m.copy();mapping=pd.read_csv(ROOT/'results/JGI_mapping/all_original_protein_mapping.csv',dtype={'JGI_gene':str})
    evidence=pd.read_csv(ROOT/'results/locus_evidence/alignment_evidence.csv')
    def gene(locus):
        sub=mapping[mapping.public_locus_tags.eq(locus)]
        assert len(sub)==1 and bool(sub.iloc[0].original_sequence_verified_in_original_genome) and bool(sub.iloc[0].original_sequence_verified_in_public_genome)
        return sub.iloc[0].JGI_gene
    edits=[('rxn90001_c0',['HMPREF3201_00478','HMPREF3201_00479','HMPREF3201_00480'],'and',
            'Bcd catalytic component plus adjacent Etf beta and alpha. Full-length reciprocal best hits to characterized M. elsdenii Bcd/Etf; remove unsupported Etf-only OR rule.'),
           ('rxn90053_c0',['HMPREF3201_00358','HMPREF3201_00359','HMPREF3201_00360','HMPREF3201_00361','HMPREF3201_00362','HMPREF3201_00363'],'and',
            'Complete contiguous RnfCDGEAB operon. Add missing RnfG and RnfE to the existing four-component AND rule. Ion stoichiometry/direction retained as assumptions.'),
           ('rxn05938_c0',['HMPREF3201_00281','HMPREF3201_00545','HMPREF3201_01026','HMPREF3201_02392'],'or',
            'Four full-length pyruvate-synthase homologs with >=50% identity and >=90% coverage to two curated PFOR references; inferred isozymes. No measured activity or enzyme capacity asserted.')]
    ledger=[];selected=[]
    for rid,loci,operator,reason in edits:
        rx=m.reactions.get_by_id(rid);before=rx.gene_reaction_rule;genes=[gene(l) for l in loci]
        rx.gene_reaction_rule=(' '+operator+' ').join(genes)
        rx.notes['v6_GPR_evidence']=reason
        rx.notes['v6_evidence_scope']='Sequence-supported computational annotation; not target-strain enzyme assay'
        for locus,gid in zip(loci,genes):
            data=mapping[mapping.JGI_gene.eq(gid)].iloc[0]
            g=m.genes.get_by_id(gid);g.name=locus;g.annotation.update({'ncbiprotein':data.public_protein_accessions})
            g.notes['public_locus_tag']=locus
            hits=evidence[evidence.locus_tag.eq(locus)].sort_values('bitscore',ascending=False).head(3)
            for h in hits.to_dict('records'):selected.append({'reaction':rid,'JGI_gene':gid,**h})
        ledger.append({'reaction':rid,'before_GPR':before,'after_GPR':rx.gene_reaction_rule,'loci':';'.join(loci),
                       'decision':reason,'stoichiometry_changed':False,'bounds_changed':False,
                       'mass_charge_imbalance':json.dumps(rx.check_mass_balance())})
    m.id='MJR8396C_sequence_supported_GPR_curated_v6'
    cobra.io.write_sbml_model(m,str(M_PATH))
    pd.DataFrame(ledger).to_csv(OUT/'GPR_edit_ledger.csv',index=False)
    pd.DataFrame(selected).to_csv(OUT/'GPR_supporting_alignment_evidence.csv',index=False)
    # Prove these are annotation edits only; stoichiometry/bounds identical.
    assert len(original.reactions)==len(m.reactions)
    for a,b in zip(original.reactions,m.reactions):
        assert a.id==b.id and a.bounds==b.bounds
        assert {x.id:v for x,v in a.metabolites.items()}=={x.id:v for x,v in b.metabolites.items()}
    other=[('acrylate_reductase','rxn90029_c0','ETF-only rule lacks a resolved catalytic reductase assignment. Preserve as unresolved; also test the reaction disabled. Homology to Bcd is not specific evidence of acryloyl reduction.'),
           ('lactate_racemase','','LarA/B/E candidates and two partial LarC-domain hits found at different loci; no complete cognate maturation system established. No racemase added.'),
           ('chain_elongation','','Core reverse-beta-oxidation homologs do not by themselves establish C5/C6 substrate specificity, termination or transport. No speculative caproate/valerate reactions added.'),
           ('CoA_transferase','rxn00875_c0','Existing reaction retained; acetate-dependent pFBA/FVA is a model mechanism, not kinetic specificity validation.'),
           ('lactoyl_CoA_dehydratase','rxn02123_c0','Existing balanced chemistry retained. Whole-complex specificity remains an inferred annotation.')]
    pd.DataFrame([{'pathway':a,'reaction':b,'decision':c} for a,b,c in other]).to_csv(OUT/'unresolved_functional_decisions.csv',index=False)
    record('results/functional_curation_v6/MJR_model_scope.json',{'source_sha256':sha(ROOT/'models/MJR8396C_sequence_annotated_primary.xml'),
           'source_gene_count':len(original.genes),'curated_gene_count':len(m.genes),'GPRs_changed':len(edits),
           'same_reaction_stoichiometry_and_bounds':True,'biomass_changed':False,'experimental_function_verified':False})
    return original,m

def verify(original,curated):
    rows=[];fvas=[];knockouts=[];cycle=[]
    variants={'M_deposited':original,'M_GPR_curated':curated,'B_deposited':load('B'),'B_DLDH_filtered':load(path=B_PATH)}
    for name,m in variants.items():
        for med in ['gut_like','carbon_rich','AA_rich']:
            medium(m,med)
            for ctl in ['baseline','no_external_lactate']:
                with m:
                    if ctl=='no_external_lactate':
                        for c in ['cpd00159','cpd00221']:
                            if 'EX_'+c+'_e0' in m.reactions:m.reactions.get_by_id('EX_'+c+'_e0').lower_bound=0
                    rows.append({'variant':name,'medium':med,'control':ctl,**phenotype(m)})
                    ids=['EX_'+c+'_e0' for c in ['cpd00159','cpd00221','cpd00211','cpd00141','cpd00029'] if 'EX_'+c+'_e0' in m.reactions]
                    f=flux_variability_analysis(m,reaction_list=ids,fraction_of_optimum=.95,processes=1).reset_index(names='reaction')
                    f['variant']=name;f['medium']=med;f['control']=ctl;fvas.append(f)
        with m:
            for r in m.boundary:r.bounds=(0,0)
            atp=m.reactions.rxn00062_c0;atp.lower_bound=0;m.objective=atp
            cycle.append({'variant':name,'ATP_max':m.slim_optimize(),'status':m.solver.status})
    for variant,m in [('before',original),('after',curated)]:
        medium(m)
        for rid in ['rxn90001_c0','rxn90053_c0','rxn05938_c0']:
            genes=sorted(g.id for g in curated.reactions.get_by_id(rid).genes)
            for gs in [[g] for g in genes]+[genes]:
                with m:
                    existing=[g for g in gs if g in m.genes]
                    cobra.manipulation.knock_out_model_genes(m,existing)
                    r=m.reactions.get_by_id(rid)
                    knockouts.append({'variant':variant,'reaction':rid,'deleted_genes':';'.join(gs),'existing_gene_count':len(existing),
                                      'reaction_disabled':r.bounds==(0,0),**phenotype(m)})
    # Reaction removal is a conservative uncertainty bound, not a functional edit.
    medium(curated)
    for rid in ['rxn90029_c0','rxn00500_c0','rxn90001_c0','rxn90053_c0','rxn05938_c0']:
        with curated:
            curated.reactions.get_by_id(rid).knock_out()
            rows.append({'variant':'M_GPR_curated','medium':'gut_like','control':'disable_'+rid,**phenotype(curated)})
    pd.DataFrame(rows).to_csv(OUT/'curation_phenotypes.csv',index=False)
    pd.concat(fvas).to_csv(OUT/'curation_exchange_FVA.csv',index=False)
    pd.DataFrame(knockouts).to_csv(OUT/'curation_gene_knockout_verification.csv',index=False)
    pd.DataFrame(cycle).to_csv(OUT/'curation_ATP_cycle_checks.csv',index=False)
    print('Functional curation and 22L evidence completed',len(rows),len(knockouts),flush=True)

if __name__=='__main__':
    OUT.mkdir(exist_ok=True);b22l_evidence();a,b=m_curation();verify(a,b)
