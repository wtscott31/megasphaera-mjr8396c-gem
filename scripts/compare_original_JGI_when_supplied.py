"""Resume exact sequence mapping only after original JGI inputs are supplied."""
from common import *
import argparse,collections,pandas as pd
from Bio import SeqIO
def main():
    p=argparse.ArgumentParser();p.add_argument('--genome',required=True);p.add_argument('--proteins',required=True);a=p.parse_args()
    genome=Path(a.genome).resolve();proteins=Path(a.proteins).resolve()
    if not genome.is_file() or not proteins.is_file():raise FileNotFoundError('Real original JGI FASTA/protein files are required; no mapping was inferred.')
    public=ROOT/'inputs/genomes/GCA_001546855.1';out=ROOT/'results/JGI_resume';out.mkdir(exist_ok=True)
    lookup=collections.defaultdict(list)
    for r in SeqIO.parse(public/'protein.faa','fasta'):lookup[str(r.seq).upper().rstrip('*')].append(r.id)
    rows=[]
    for r in SeqIO.parse(proteins,'fasta'):
        hits=lookup[str(r.seq).upper().rstrip('*')]
        rows.append({'original_JGI_record':r.id,'public_protein_accessions':';'.join(hits),'status':'unique_exact_protein_match' if len(hits)==1 else ('ambiguous_identical_sequences' if hits else 'no_exact_match'),'length_AA':len(r.seq)})
    pd.DataFrame(rows).to_csv(out/'exact_protein_mapping.csv',index=False)
    record('results/JGI_resume/input_provenance.json',{'genome':str(genome),'genome_sha256':sha(genome),'proteins':str(proteins),'proteins_sha256':sha(proteins),'run_utc':now(),'scope':'Exact amino-acid sequence equality only; identical paralogs remain ambiguous. No nonexact GPR mapping is promoted.'})
    for q,r,label in [(genome,public/'genomic.fna','JGI_vs_public'),(public/'genomic.fna',genome,'public_vs_JGI')]:run(['skani','dist','--slow','-t','4','-q',q,'-r',r,'-o',out/(label+'.tsv')],'JGI_resume/'+label)
if __name__=='__main__':main()
