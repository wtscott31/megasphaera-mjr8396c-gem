"""Sequence statistics, checksum validation, exact GTDB lookup, reciprocal ANI."""
from common import *
import re,gzip,collections,io
from Bio import SeqIO
import pandas as pd

def n50(lengths):
    lengths=sorted(lengths,reverse=True);t=sum(lengths)/2;s=0
    for v in lengths:
        s+=v
        if s>=t:return v
    return 0
def stats():
    rows=[];features=[];valid=[]
    for folder in sorted((ROOT/'inputs/genomes').iterdir()):
        f=folder/'genomic.fna'
        if not f.exists():continue
        rec=list(SeqIO.parse(f,'fasta'));seq=''.join(str(r.seq).upper() for r in rec);c=collections.Counter(seq)
        contigs=[x for r in rec for x in re.split('N+',str(r.seq).upper()) if x]
        row={'accession':folder.name,'fasta_sha256':sha(f),'total_bp':len(seq),'canonical_bp':sum(c[x] for x in 'ACGT'),'ambiguous_bp':sum(c[x] for x in c if x not in 'ACGT'),'fasta_records_scaffolds':len(rec),'scaffold_N50_bp':n50([len(r) for r in rec]),'N_split_segments':len(contigs),'N_split_contig_N50_bp':n50([len(x) for x in contigs]),'GC_percent_canonical':100*(c['G']+c['C'])/sum(c[x] for x in 'ACGT'),**{x:c[x] for x in 'ACGTN'},'fasta_ids':';'.join(r.id for r in rec)}
        gb=folder/'genomic.gbff'
        if gb.exists():
            gr=list(SeqIO.parse(gb,'genbank'));row['genbank_records']=len(gr);row['genbank_organisms']=';'.join(sorted({r.annotations.get('organism','') for r in gr}))
            counts=collections.Counter(x.type for r in gr for x in r.features);row.update({f'features_{k}':v for k,v in counts.items()})
            for r in gr:
                for feat in r.features:
                    if feat.type!='CDS':continue
                    q=feat.qualifiers
                    features.append({'assembly':folder.name,'scaffold':r.id,'start_1based':int(feat.location.start)+1,'end_1based_inclusive':int(feat.location.end),'strand':feat.location.strand,'locus_tag':';'.join(q.get('locus_tag',[])),'protein_id':';'.join(q.get('protein_id',[])),'gene':';'.join(q.get('gene',[])),'product':';'.join(q.get('product',[])),'pseudo':'pseudo' in q or 'pseudogene' in q,'translation_length':len(q.get('translation',[''])[0]),'EC':';'.join(q.get('EC_number',[]))})
            # Scaffold FASTA must agree with annotated GenBank sequences.
            row['fasta_genbank_sequence_multiset_equal']=sorted(str(r.seq).upper() for r in rec)==sorted(str(r.seq).upper() for r in gr)
        faa=folder/'protein.faa'
        if faa.exists():row['protein_fasta_records']=sum(1 for _ in SeqIO.parse(faa,'fasta'))
        md=folder/'md5checksums.txt'
        if md.exists():
            expected={Path(line.split()[1]).name:line.split()[0] for line in md.read_text().splitlines() if len(line.split())>=2}
            meta=json.loads((ROOT/f'provenance/assemblies/{folder.name}.json').read_text())['reports'][0]
            prefix=folder.name+'_'+meta['assembly_info']['assembly_name']+'_'
            for p in folder.glob('*.gz'):
                wanted=expected.get(prefix+p.name)
                actual=hashlib.md5(p.read_bytes()).hexdigest();valid.append({'assembly':folder.name,'file':p.name,'expected_md5':wanted,'actual_md5':actual,'matches':actual==wanted})
        rows.append(row)
    pd.DataFrame(rows).to_csv(ROOT/'results/genome_statistics.csv',index=False)
    pd.DataFrame(features).to_csv(ROOT/'results/genbank_CDS_locus_ledger.csv',index=False)
    pd.DataFrame(valid).to_csv(ROOT/'results/NCBI_md5_validation.csv',index=False)
    a=ROOT/'inputs/genomes/GCA_001546855.1/genomic.fna';b=ROOT/'inputs/genomes/GCF_001546855.1/genomic.fna'
    record('results/target_assembly_validation.json',{'genbank_refseq_same_sequence_multiset':sorted(str(x.seq).upper() for x in SeqIO.parse(a,'fasta'))==sorted(str(x.seq).upper() for x in SeqIO.parse(b,'fasta')),'local_JGI_comparison':'Pending original JGI sequence input; numeric model IDs do not establish equivalence.','assembly_limitation':'Draft scaffolds contain assembly Ns; N-split contigs are inferred from sequence runs, not read-supported closure.'})
    return rows
def taxonomy():
    p=ROOT/'inputs/gtdb/bac120_taxonomy_r232.tsv.gz';rows=[]
    with gzip.open(p,'rt') as f:
        for line in f:
            if 'Megasphaera' in line:rows.append(line.rstrip().split('\t'))
    pd.DataFrame(rows,columns=['accession','taxonomy']).to_csv(ROOT/'results/gtdb_Megasphaera_taxonomy_r232.csv',index=False)
    hit=[r for r in rows if '001546855.1' in r[0]]
    record('results/GTDB_exact_assembly_lookup.json',{'method':'exact assembly lookup; NOT a new GTDB-Tk classification','release':'GTDB R232 / 232.0','table_sha256':sha(p),'NCBI_GenBank_RefSeq_pair_metadata':'provenance/assemblies/GCA_001546855.1.json','result':hit,'GTDB_Tk_status':'not_run','reason':'No compatible GTDB-Tk reference database or supported Linux execution environment supplied. Host is macOS arm64 with 64 GiB RAM. Only the small taxonomy/tree resources were downloaded.'})
def ani():
    fastas=sorted((ROOT/'inputs/genomes').glob('*/genomic.fna'))
    # Target GenBank versus every available assembly, both query directions.
    target=ROOT/'inputs/genomes/GCA_001546855.1/genomic.fna'; rows=[]
    out=ROOT/'results/ani';out.mkdir(exist_ok=True)
    for ref in fastas:
        if ref==target:continue
        for q,r in [(target,ref),(ref,target)]:
            label=q.parent.name+'__'+r.parent.name
            dest=out/(label+'.skani.tsv')
            result=run(['skani','dist','--slow','--min-af','15','--ci','-t','4','-q',q,'-r',r,'-o',dest],'ani/skani_'+label,180)
            if dest.exists() and len(dest.read_text().splitlines())>1:
                raw=pd.read_csv(dest,sep='\t').to_dict('records')[0];rows.append({'query':q.parent.name,'reference':r.parent.name,'tool':'skani 0.3.1','status':'estimate_returned',**raw})
            else:rows.append({'query':q.parent.name,'reference':r.parent.name,'tool':'skani 0.3.1','status':'no_estimate','meaning':'No estimate returned at the recorded screening/alignment settings; identity is not zero.'})
            dest=out/(label+'.fastANI.tsv');run(['fastANI','-q',q,'-r',r,'-t','4','-o',dest],'ani/fastANI_'+label,180)
            if dest.exists() and dest.stat().st_size:
                parts=dest.read_text().strip().split('\t');rows.append({'query':q.parent.name,'reference':r.parent.name,'tool':'FastANI 1.34','status':'estimate_returned','ANI':float(parts[2]),'mapped_fragments':int(parts[3]),'total_query_fragments':int(parts[4]),'query_aligned_fraction':int(parts[3])/int(parts[4])})
            else:rows.append({'query':q.parent.name,'reference':r.parent.name,'tool':'FastANI 1.34','status':'no_estimate','meaning':'No estimate; do not report 0% ANI or infer species from absence.'})
        pd.DataFrame(rows).to_csv(ROOT/'results/reciprocal_ANI.csv',index=False)
if __name__=='__main__':
    stats();taxonomy();ani()
