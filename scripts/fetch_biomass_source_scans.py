"""Re-fetch official primary-source scans into a local research cache.

Full non-open-access articles are excluded from the distributed release.
Extracted numerical facts and bibliographic/retrieval provenance are included.
"""
from common import *
import re

def main():
    for pmcid in ['239026','216544','243251']:
        p=fetch('https://pmc.ncbi.nlm.nih.gov/articles/PMC'+pmcid+'/?report=xml','inputs/biomass_literature/PMC'+pmcid+'.html',30)
        if p is None:continue
        images=list(dict.fromkeys(re.findall(r'(?:src|href)="(https://cdn\.ncbi\.nlm\.nih\.gov/pmc/blobs/[^" ]+\.png)"',p.read_text())))
        if not images:raise RuntimeError('No primary-source page scans resolved for PMC'+pmcid)
        for i,url in enumerate(images):fetch(url,f'inputs/biomass_literature/PMC{pmcid}_page{i+1}.png',30)
if __name__=='__main__':main()
