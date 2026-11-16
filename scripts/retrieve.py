"""Download versioned public inputs; every attempt is recorded, including failures."""
from common import *
import gzip, argparse

def assembly(acc, full=False):
    metadata=ROOT/f'provenance/assemblies/{acc}.json'
    if not metadata.exists():
        if not fetch(f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/dataset_report',str(metadata.relative_to(ROOT))):return
    data=json.loads(metadata.read_text()); reports=data.get('reports',[])
    if not reports: return
    rep=reports[0]; name=rep['assembly_info']['assembly_name']; stem=acc+'_'+name
    prefix,number=acc.split('_'); number=number.split('.')[0]
    url=f'https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{number[:3]}/{number[3:6]}/{number[6:9]}/{stem}/{stem}'
    extensions=['genomic.fna.gz','assembly_report.txt','md5checksums.txt']
    if full:extensions+=['genomic.gbff.gz','protein.faa.gz','genomic.gff.gz']
    for ext in extensions:
        suffix=ext if ext=='md5checksums.txt' else '_'+ext
        link=(url.rsplit('/',1)[0]+'/md5checksums.txt') if ext=='md5checksums.txt' else url+suffix
        dest=f'inputs/genomes/{acc}/{ext}'
        p=ROOT/dest
        if not p.exists():p=fetch(link,dest,90)
        if p and ext.endswith('.gz'):
            try:
                raw=gzip.decompress(p.read_bytes());p.with_suffix('').write_bytes(raw)
            except Exception as e:record(f'logs/{acc}_{ext}_error.json',{'error':repr(e)})
    return rep

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('accessions',nargs='*');p.add_argument('--full',action='store_true');p.add_argument('--resources',action='store_true');a=p.parse_args()
    for acc in a.accessions:assembly(acc,a.full)
    if a.resources:
        fetch('https://data.gtdb.ecogenomic.org/releases/release232/232.0/','provenance/gtdb_release232_index.html')
        fetch('https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_taxonomy_r232.tsv.gz','inputs/gtdb/bac120_taxonomy_r232.tsv.gz',120)
        fetch('https://pmc.ncbi.nlm.nih.gov/articles/PMC7232489/?report=xml','provenance/iME375_article.html')
        fetch('https://www.mdpi.com/2076-2607/8/4/539/s1','provenance/iME375_supplement_index.html')
        fetch('https://www.vmh.life/files/reconstructions/AGORA2/version2.01/sbml_files/individual_reconstructions/','provenance/agora2_index.html')
        fetch('https://www.vmh.life/files/reconstructions/AGORA2/version2.01/README.txt','provenance/AGORA2_README.txt')
        fetch('https://raw.githubusercontent.com/VirtualMetabolicHuman/AGORA2/master/LICENSE','provenance/AGORA2_LICENSE.txt')
