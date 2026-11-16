"""NEW uniform gapseq reconstructions; never substitutions for published GEMs.

Run using the read-only existing gapseq-clean runtime. All outputs and copied
gapseq source/reference data are isolated under this release directory.
"""
from common import *
import os, shutil, concurrent.futures, argparse
import pandas as pd
from Bio import SeqIO

STRAINS={'MJR8396C':('GCA_001546855.1','neg'),
         'NM10':('GCA_000417505.1','neg'),
         'BL7':('GCA_000417525.1','neg'),
         'B22L':('GCF_000737885.1','pos')}
ENV_PATH=os.environ.get('PATH','')

def prepare():
    out=ROOT/'results/matched_reconstruction_v6';out.mkdir(exist_ok=True)
    media=pd.read_csv(BASE/'media/MJR8396C_gut_like_medium.tsv',sep='\t')
    media[['compounds','name','maxFlux']].to_csv(out/'shared_gapfill_medium.tsv',sep='\t',index=False)
    records=[]
    for name,(accession,gram) in STRAINS.items():
        d=out/name;d.mkdir(exist_ok=True)
        source=ROOT/'inputs/genomes'/accession/'genomic.fna'
        if not source.exists() and name=='B22L':
            gb=source.with_name('genomic.gbff')
            SeqIO.write(SeqIO.parse(gb,'genbank'),source,'fasta')
        target=d/(name+'.fna');shutil.copy2(source,target)
        records.append({'name':name,'accession':accession,'biomass_template':gram,
                        'source':str(source.relative_to(ROOT)),'source_sha256':sha(source),
                        'input':str(target.relative_to(ROOT)),
                        'annotation':'uniform native gapseq pyrodigal translation; NOT published NM10/BL7 models'})
    record('provenance/matched_reconstruction_v6_inputs.json',records)

def reconstruct(name):
    accession,gram=STRAINS[name]
    d=ROOT/'results/matched_reconstruction_v6'/name
    gp=os.environ.get('GAPSEQ_BIN') or shutil.which('gapseq');fna=d/(name+'.fna')
    if not gp: raise FileNotFoundError('Install pinned gapseq 2.1.0 and set GAPSEQ_BIN or PATH; see docs/REPRODUCIBILITY.md')
    steps=[('find',[str(gp),'find','-p','all','-A','diamond','-K','4','-t','Bacteria','-M','nucl','-b','200','-c','75','-O','-f',str(d),str(fna)],d/(name+'-all-Reactions.tbl')),
           ('transport',[str(gp),'find-transport','-A','diamond','-K','4','-M','nucl','-b','50','-c','75','-f',str(d),str(fna)],d/(name+'-Transporter.tbl')),
           ('draft',[str(gp),'draft','-r',str(d/(name+'-all-Reactions.tbl')),'-t',str(d/(name+'-Transporter.tbl')),'-p',str(d/(name+'-all-Pathways.tbl')),'-b',gram,'-n',name,'-f',str(d)],d/(name+'-draft.RDS')),
           ('fill',[str(gp),'fill','-m',str(d/(name+'-draft.RDS')),'-n',str(d.parent/'shared_gapfill_medium.tsv'),'-z','glpk','-k','0.01','-f',str(d)],d/(name+'.RDS'))]
    history=[]
    env=os.environ.copy();env.update(PATH=ENV_PATH,R_LIBS_USER='/tmp/MJR8396C_v6_empty_R',R_ENVIRON_USER='/dev/null',R_PROFILE_USER='/dev/null',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    for stage,cmd,expected in steps:
        status=d/(stage+'_status.json')
        if status.exists() and json.loads(status.read_text()).get('returncode')==0 and expected.exists():
            history.append(json.loads(status.read_text()));continue
        log=d/(stage+'.log');start=time.monotonic()
        row={'strain':name,'accession':accession,'stage':stage,'command':cmd,
             'cwd':str(d),'started_utc':now(),'log':str(log.relative_to(ROOT)),
             'runtime':'gapseq-clean read-only existing conda environment',
             'gapseq_version':'2.1.0','sequence_database':'Bacteria 1.5 Zenodo 20446806',
             'new_reconstruction':True}
        print(name,stage,'started',now(),flush=True)
        with log.open('w') as f:
            f.write(json.dumps(row)+'\n');f.flush()
            try:
                p=subprocess.run(cmd,cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=14400)
                row['returncode']=p.returncode
            except Exception as e:row['error']=repr(e)
        row.update(seconds=time.monotonic()-start,finished_utc=now(),expected_output_exists=expected.exists())
        status.write_text(json.dumps(row,indent=2));history.append(row)
        print(name,stage,'finished',row.get('returncode',row.get('error')),round(row['seconds'],1),flush=True)
        if row.get('returncode')!=0 or not expected.exists():break
    (d/'execution_history.json').write_text(json.dumps(history,indent=2))
    return history

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--strains',nargs='+',default=list(STRAINS));args=parser.parse_args()
    prepare()
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        results=list(pool.map(reconstruct,args.strains))
    record('results/matched_reconstruction_v6/execution_summary.json',results)
