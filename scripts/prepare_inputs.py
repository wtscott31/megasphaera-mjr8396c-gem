"""Rebuild derived input adapters from preserved raw files, without downloads."""
from common import *
import zipfile,collections
from Bio import SeqIO
from lxml import etree
def main():
    z=ROOT/'inputs/comparators/iME375_mdpi_supplement.zip'
    with zipfile.ZipFile(z) as f:
        for n in f.namelist():
            if Path(n).is_absolute() or '..' in Path(n).parts:raise ValueError('Unsafe ZIP path')
        f.extractall(ROOT/'inputs/comparators/iME375')
    refs=ROOT/'inputs/reference_proteins';entries={};sequences={}
    for stem in ['reviewed_energy_lactate','chain_elongation_additional']:
        for x in json.loads((refs/(stem+'.json')).read_text())['results']:entries[x['primaryAccession']]=x
        for x in SeqIO.parse(refs/(stem+'.faa'),'fasta'):sequences[x.id]=x
    record('inputs/reference_proteins/combined_reference_panel.json',{'results':list(entries.values())})
    SeqIO.write(sequences.values(),refs/'combined_reference_panel.faa','fasta')
    rows=[]
    for p in (ROOT/'inputs/comparators/AGORA2_fixed').glob('*.xml'):
        data=p.read_bytes()
        try:etree.fromstring(data);continue
        except Exception as e:error=str(e)
        text=data.decode('utf-8',errors='surrogateescape')
        bad=[{'character_offset':i,'invalid_byte_hex':hex(ord(c)-0xdc00),'mapped_latin1_character':bytes([ord(c)-0xdc00]).decode('latin1')} for i,c in enumerate(text) if 0xdc80<=ord(c)<=0xdcff]
        if not bad:raise ValueError(error)
        fixed=''.join(bytes([ord(c)-0xdc00]).decode('latin1') if 0xdc80<=ord(c)<=0xdcff else c for c in text).encode('utf-8');etree.fromstring(fixed)
        dest=ROOT/'inputs/comparators/encoding_repaired'/p.name;dest.parent.mkdir(exist_ok=True);dest.write_bytes(fixed)
        rows.append({'original':str(p.relative_to(ROOT)),'derived':str(dest.relative_to(ROOT)),'original_sha256':sha(p),'derived_sha256':sha(dest),'original_error':error,'only_change':'Invalid non-UTF8 annotation bytes mapped via Latin1 to valid UTF8. No chemistry, identifier, bound or GPR changed.','replacements':bad})
    record('provenance/AGORA2_annotation_encoding_repairs.json',rows)
if __name__=='__main__':main()
