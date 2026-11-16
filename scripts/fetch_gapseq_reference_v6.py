"""Restore the exact excluded public gapseq Bacteria reference archive."""
from common import *
import tarfile

URL='https://zenodo.org/api/records/20446806/files/Bacteria.tar.gz/content'
MD5='09db671cb466f6bab5edac698bd92d81'
SHA256='ec669d568a2e5459bfe6619865a3cd2612fb9edc3a9688ad45fdf619d2e1afeb'
def unpack(archive,dest):
    with tarfile.open(archive) as t:
        members=t.getmembers()
        for m in members:
            p=Path(m.name)
            if p.is_absolute() or '..' in p.parts or m.issym() or m.islnk() or not (m.isfile() or m.isdir()):raise ValueError('Unsafe archive member '+m.name)
        t.extractall(dest,members=members)
def main():
    p=ROOT/'inputs/Bacteria_gapseq_DB_1.5_20446806.tar.gz'
    if not p.exists():
        with (ROOT/'logs/v6_reference_restore.log').open('w') as log:
            subprocess.run(['/usr/bin/curl','--http1.1','-L','--fail','--show-error','--max-time','600',URL,'-o',str(p)],stdout=log,stderr=subprocess.STDOUT,check=True)
    assert hashlib.md5(p.read_bytes()).hexdigest()==MD5 and sha(p)==SHA256
    dest=ROOT/'tools/gapseq/dat/seq';dest.mkdir(parents=True,exist_ok=True);unpack(p,dest)
    for sub in ['rev','unrev','rxn']:unpack(dest/'Bacteria'/sub/'sequences.tar.gz',dest/'Bacteria'/sub)
    (dest/'Bacteria/version_seqDB.json').write_text(json.dumps({'zenodoID':20446806,'version':'1.5','date':'2026-05-29'}))
    record('provenance/v6_reference_restored.json',{'restored_utc':now(),'url':URL,'sha256':SHA256,'md5':MD5,'verified':True})
    print('Pinned Bacteria1.5 references restored and verified.')
if __name__=='__main__':main()
