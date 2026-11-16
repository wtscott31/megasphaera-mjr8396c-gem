"""Sequence and genomic-coordinate mapping of the actual original IMG proteins."""
from native_analysis import *
from Bio.Data import CodonTable
from Bio.Seq import Seq

def six_frames(records):
    rows=[]
    for rec in records:
        for strand,seq in [(1,rec.seq),(-1,rec.seq.reverse_complement())]:
            for frame in range(3):
                size=((len(seq)-frame)//3)*3
                rows.append((rec.id,strand,frame,len(rec.seq),str(seq[frame:frame+size].translate(table=11)),str(seq)))
    return rows

def positions(protein,frames):
    found=set();protein=protein.rstrip('*');starts=set(CodonTable.unambiguous_dna_by_id[11].start_codons)
    for rid,strand,frame,n,aa,dna in frames:
        # Bacterial alternative initiation is verified from the actual codon.
        query=protein[1:] if protein.startswith('M') else protein
        skip=1 if protein.startswith('M') else 0
        start=aa.find(query)
        while start>=0:
            idx=start-skip
            if idx>=0:
                pos=frame+idx*3;end=pos+len(protein)*3
                if skip==0 or aa[idx]=='M' or dna[pos:pos+3] in starts:
                    lo,hi=(pos,end) if strand==1 else (n-end,n-pos)
                    found.add((rid,lo,hi,strand,'exact_translation' if aa[idx:idx+len(protein)]==protein else 'verified_alternative_start_translation'))
            start=aa.find(query,start+1)
    return sorted(found)

def main():
    out=ROOT/'results/JGI_mapping';out.mkdir(exist_ok=True)
    jpath=ROOT/'inputs/original_JGI_gapseq/Megasphaera_genome_assembly'
    jr=list(SeqIO.parse(jpath/'2775506928.fna','fasta'));pr=list(SeqIO.parse(ROOT/'inputs/genomes/GCA_001546855.1/genomic.fna','fasta'))
    proteins=list(SeqIO.parse(jpath/'2775506928_genes.faa','fasta'))
    jf=six_frames(jr);pf=six_frames(pr)
    pubs=list(SeqIO.parse(ROOT/'inputs/genomes/GCA_001546855.1/protein.faa','fasta'))
    exact=collections.defaultdict(list)
    for r in pubs:exact[str(r.seq).rstrip('*')].append(r.id)
    features=[];byprot={}
    for r in SeqIO.parse(ROOT/'inputs/genomes/GCA_001546855.1/genomic.gbff','genbank'):
        for ft in r.features:
            if ft.type!='CDS' or 'protein_id' not in ft.qualifiers:continue
            pid=ft.qualifiers['protein_id'][0];lo=int(ft.location.start);hi=int(ft.location.end);strand=ft.location.strand
            coding=ft.extract(r.seq)
            stop=str(coding[-3:]).upper() in ['TAA','TAG','TGA']
            if stop:
                if strand==1:hi-=3
                else:lo+=3
            row={'scaffold':r.id,'start':lo,'end':hi,'strand':strand,'protein':pid,'locus':ft.qualifiers.get('locus_tag',[''])[0],'product':ft.qualifiers.get('product',[''])[0]}
            features.append(row);byprot[pid]=row
    rows=[];genomic=[]
    for i,r in enumerate(proteins):
        seq=str(r.seq).rstrip('*');jp=positions(seq,jf);pp=positions(seq,pf)
        matches=exact[seq]
        candidates=[]
        if not matches:
            for rid,lo,hi,strand,kind in pp:
                for ft in features:
                    if ft['scaffold']!=rid or ft['strand']!=strand:continue
                    overlap=max(0,min(hi,ft['end'])-max(lo,ft['start']))
                    # The same exact translated genomic ORF, same reading frame
                    # and same protein C terminus establishes a boundary variant.
                    same_end=(hi==ft['end']) if strand==1 else (lo==ft['start'])
                    if same_end and overlap>=min(hi-lo,ft['end']-ft['start'])*.95 and (lo-ft['start'])%3==0:
                        candidates.append(ft['protein'])
            matches=sorted(set(candidates))
            status='unique_coordinate_verified_annotation_boundary_variant' if len(matches)==1 else ('ambiguous_coordinate_overlap' if matches else 'no_equivalent_public_CDS')
        else:status='unique_exact_protein_match' if len(matches)==1 else 'ambiguous_identical_sequences'
        rows.append({'JGI_gene':r.id,'original_header':r.description,'length_AA':len(seq),'AA_sha256':hashlib.sha256(seq.encode()).hexdigest(),'original_genomic_locations':len(jp),'public_genomic_locations':len(pp),'public_protein_accessions':';'.join(matches),'public_locus_tags':';'.join(byprot[x]['locus'] for x in matches),'mapping_status':status,'original_sequence_verified_in_original_genome':bool(jp),'original_sequence_verified_in_public_genome':bool(pp)})
        for origin,locs in [('original_JGI',jp),('public_GCA_001546855.1',pp)]:
            for rid,lo,hi,strand,kind in locs:genomic.append({'JGI_gene':r.id,'assembly':origin,'scaffold':rid,'start_0based':lo,'end_exclusive':hi,'strand':strand,'translation_evidence':kind})
        if i%250==0:print('Mapped',i,'of',len(proteins),flush=True)
    mapping=pd.DataFrame(rows);mapping.to_csv(out/'all_original_protein_mapping.csv',index=False)
    pd.DataFrame(genomic).to_csv(out/'all_genomic_coordinates.csv',index=False)
    m=load();gids={g.id for g in m.genes};used=mapping[mapping.JGI_gene.isin(gids)].copy();used.to_csv(out/'deposited_model_477_gene_mapping.csv',index=False)
    missing=gids-set(mapping.JGI_gene);assert not missing,missing
    reaction=[]
    byid=used.set_index('JGI_gene').to_dict('index')
    for rxn in m.reactions:
        for g in rxn.genes:
            reaction.append({'reaction':rxn.id,'name':rxn.name,'current_GPR':rxn.gene_reaction_rule,'JGI_gene':g.id,**byid[g.id],'decision':'Retain original GPR; sequence identity/coordinate mapping establishes locus, not enzyme specificity'})
    pd.DataFrame(reaction).to_csv(out/'reaction_GPR_public_locus_ledger.csv',index=False)
    # Determine which original contigs are literally present in public bytes.
    contigs=[]
    for r in jr:
        seq=str(r.seq).upper();hits=[]
        for p in pr:
            for strand,s in [(1,seq),(-1,str(Seq(seq).reverse_complement()))]:
                start=str(p.seq).upper().find(s)
                if start>=0:hits.append((p.id,start,start+len(s),strand))
        contigs.append({'JGI_scaffold':r.id,'length':len(seq),'exact_public_locations':json.dumps(hits),'exact_match':bool(hits)})
    pd.DataFrame(contigs).to_csv(out/'exact_contig_comparison.csv',index=False)
    orig=collections.Counter(''.join(str(r.seq).upper() for r in jr));pub=collections.Counter(''.join(str(r.seq).upper() for r in pr))
    summary={'original_scaffolds':len(jr),'original_total_bp':sum(map(lambda r:len(r.seq),jr)),'public_scaffolds':len(pr),'public_total_bp':sum(map(lambda r:len(r.seq),pr)),'identical_assembly':False,'original_ambiguous_bases':sum(v for k,v in orig.items() if k not in 'ACGT'),'original_GC_percent_canonical':100*(orig['G']+orig['C'])/sum(orig[x] for x in 'ACGT'),'original_protein_count':len(proteins),'all_mapping_status_counts':mapping.mapping_status.value_counts().to_dict(),'model_gene_mapping_status_counts':used.mapping_status.value_counts().to_dict(),'model_genes_verified_in_original_genome':int(used.original_sequence_verified_in_original_genome.sum()),'model_genes_verified_in_public_genome':int(used.original_sequence_verified_in_public_genome.sum()),'model_genes_total':len(gids),'original_contigs_exactly_present_in_public':sum(x['exact_match'] for x in contigs),'header_typo':'Original IMG files say MJR9396C; sequence comparison, not the label, supports same strain lineage','interpretation':'100% skani ANI over matched sequence is not identity of these differently fragmented/filtered assemblies. Public annotation and original gene calls remain separate.'}
    record('results/JGI_mapping/summary.json',summary);print(summary,flush=True)

if __name__=='__main__':main()
