from native_analysis import *
import argparse
def main(name,solver='hybrid'):
    m=load(path=ROOT/'results/matched_reconstruction_v6'/name/(name+'.xml'))
    m.solver=solver;m.solver.configuration.timeout=15
    out=ROOT/'results/matched_monocultures_v6'/name
    rows=[]
    for med in ['gut_like','carbon_rich','AA_rich']:
        medium(m,med)
        ids=['EX_'+c+'_e0' for c in ['cpd00211','cpd00141','cpd00029','cpd00159','cpd00221','cpd00239','cpd00324'] if 'EX_'+c+'_e0' in m.reactions]
        start=time.monotonic()
        try:
            f=flux_variability_analysis(m,reaction_list=ids,fraction_of_optimum=.95,loopless='fastSNP',processes=1).reset_index(names='reaction')
            for x in f.to_dict('records'):rows.append({'model':name,'medium':med,'method':'native COBRApy fastSNP','growth_fraction':.95,'status':'optimal','seconds':time.monotonic()-start,**x})
        except Exception as e:rows.append({'model':name,'medium':med,'method':'native COBRApy fastSNP','status':'failed','error':repr(e),'seconds':time.monotonic()-start})
        pd.DataFrame(rows).to_csv(out/('loopless_FVA95_'+solver+'.csv'),index=False)
        print(name,med,'native fastSNP completed',time.monotonic()-start,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--strain',required=True);p.add_argument('--solver',default='hybrid');a=p.parse_args();main(a.strain,a.solver)
