"""Extract selected original AGORA2 ZIP members using HTTP byte ranges."""
from common import *
import io,zipfile,re,urllib.parse,struct,zlib,copy
class RemoteZip(io.RawIOBase):
    def __init__(self,url,size,label):self.url=url;self.size=size;self.pos=0;self.label=label;self.calls=0
    def seekable(self):return True
    def readable(self):return True
    def tell(self):return self.pos
    def seek(self,offset,whence=0):
        self.pos=offset if whence==0 else (self.pos+offset if whence==1 else self.size+offset);return self.pos
    def read(self,n=-1):
        n=self.size-self.pos if n<0 else min(n,self.size-self.pos)
        if n<=0:return b''
        if n>20_000_000:raise ValueError('Refuse range over 20 MB; no bulk archive download.')
        start=self.pos;end=self.pos+n-1;self.calls+=1
        p=ROOT/f'inputs/comparators/archive_ranges/{self.label}_{start}_{end}.bin';p.parent.mkdir(parents=True,exist_ok=True)
        cmd=['/usr/bin/curl','--http1.1','-L','--fail','--silent','--show-error','--max-time','90','--max-filesize','20000000','--range',f'{start}-{end}','-D',str(p)+'.headers',self.url,'-o',str(p)]
        if not p.exists():
            r=run(cmd,'remote_zip/'+self.label+'_'+str(self.calls),100)
            if r.get('returncode')!=0:raise IOError(r)
        b=p.read_bytes();headers=Path(str(p)+'.headers').read_text()
        if len(b)!=n or f'bytes {start}-{end}/{self.size}' not in headers:raise IOError('Server did not return the requested verified byte range.')
        self.pos+=len(b);return b
def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--fixed',action='store_true');args=parser.parse_args()
    path='AGORA2/version2.01/sbml_files_fixed/zipped/AGORA2_models/AGORA2_SBML.zip' if args.fixed else 'AGORA2/version2.01/sbml_files/zipped/AGORA2_annotatedSBML_all.zip'
    label='fixed_all' if args.fixed else 'original_all';folder='AGORA2_fixed' if args.fixed else 'AGORA2_original'
    url='https://www.vmh.life/api/reconstructionFiles?download='+urllib.parse.quote(path)
    raw=RemoteZip(url,2170499647 if args.fixed else 6202657544,label)
    with zipfile.ZipFile(raw) as z:
        names=[i for i in z.infolist() if 'Megasphaera' in i.filename and i.filename.endswith('.xml')]
        record('provenance/'+folder+'_ZIP_selected_members.json',[{'filename':i.filename,'CRC':i.CRC,'file_size':i.file_size,'compress_size':i.compress_size,'header_offset':i.header_offset} for i in names])
        for i in names:
            # Python 3.11 lacks the archive's Deflate64 compression. Preserve
            # each original compressed member in a standalone ZIP and decode
            # with macOS libarchive/bsdtar; validate original CRC and byte count.
            raw.seek(i.header_offset);header=raw.read(30);fields=struct.unpack('<4s5H3I2H',header)
            payload=header+raw.read(fields[-2]+fields[-1]+i.compress_size)
            single=ROOT/'inputs/comparators/archive_ranges'/(label+'_'+Path(i.filename).stem+'.zip')
            buffer=io.BytesIO()
            with zipfile.ZipFile(buffer,'w') as zz:
                zz.fp.write(payload);info=copy.copy(i);info.header_offset=0
                zz.filelist=[info];zz.start_dir=zz.fp.tell();zz._didModify=True
            single.write_bytes(buffer.getvalue())
            result=subprocess.run(['/usr/bin/bsdtar','-xOf',str(single),i.filename],capture_output=True)
            if result.returncode:
                if i.compress_type!=9:raise RuntimeError(result.stderr.decode(errors='replace'))
                import inflate64
                data=inflate64.Inflater().inflate(payload[30+fields[-2]+fields[-1]:])
            else:data=result.stdout
            if len(data)!=i.file_size or zlib.crc32(data)!=i.CRC:raise IOError('Original ZIP member CRC/length mismatch')
            dest=ROOT/'inputs/comparators'/folder/Path(i.filename).name;dest.parent.mkdir(exist_ok=True);dest.write_bytes(data)
            print(dest.name,len(data),flush=True)
    record('provenance/'+folder+'_remote_archive_provenance.json',{'url':url,'archive_bytes':raw.size,'method':'Verified HTTP byte-range reads; original ZIP member CRC and length verified explicitly. Whole archive was not downloaded or SHA256-hashed; selected raw ranges and extracted-model SHA256s are preserved.','members':{p.name:sha(p) for p in (ROOT/'inputs/comparators'/folder).glob('*.xml')}})
if __name__=='__main__':main()
