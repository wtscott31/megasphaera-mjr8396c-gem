"""Export exact original notebook cells and outputs; do not execute heavy MILPs."""
from common import *
import re,pandas as pd
def main():
    folder=ROOT/'inputs/original_project';rows=[];seeds=[]
    for p in sorted(folder.glob('*.ipynb')):
        n=json.loads(p.read_text());cells=[]
        for i,c in enumerate(n['cells']):
            text=''.join(c.get('source',[]))
            if re.search('SDModule|compute_strain_designs|tilt|ko_cost|solver|objective|candidate|time_limit|Community|tradeoff',text,re.I):
                cells.append({'cell_index_zero_based':i,'execution_count':c.get('execution_count'),'source':text,'outputs':c.get('outputs',[])})
            for j,o in enumerate(c.get('outputs',[])):
                out=''.join(o.get('text',[]))
                seeds.extend({'notebook':p.name,'cell_index':i,'seed':int(seed)} for seed in re.findall('random seed (\\d+)',out))
        record('provenance/'+p.stem+'_verbatim_source_and_outputs.json',{'source_sha256':sha(p),'cells':cells})
    record('results/recovered_random_seeds.json',seeds)
    p=folder/'outputs/tables/strain_design_priority_reaction_candidates.csv';cost=pd.read_csv(p);cost['recovered_KO_cost']=1.
    cost.to_csv(ROOT/'results/recovered_original_78_candidate_costs.csv',index=False)
    assert len(cost)==78
    # Full notebook and earlier detailed provenance JSON remain immutable
    # evidence for the recovery narrative and missing version/proof metadata.
if __name__=='__main__':main()
