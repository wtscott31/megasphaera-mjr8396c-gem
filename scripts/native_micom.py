"""Native MICOM controls with open-source HiGHS/OSQP via optlang hybrid."""
from native_analysis import *
from micom import Community
import traceback

def main():
    out=ROOT/'results/native_micom';out.mkdir(exist_ok=True)
    settings={'engine':'native MICOM 0.39.1','solver':'optlang hybrid / HiGHS LP + OSQP QP','abundances':[.1,.25,.5,.75,.9],'fractions':[.5,.8,.95,1.],'max_member_import':100,'min_growth':0,'atol':1e-3,'rtol':1e-3,'native_source_modified':False,'solver_time_limit_seconds':5,'notes':'Native QP/pFBA numerical sensitivity run. Growth preservation tolerance 0.001 absolute/relative is explicit, after tighter native attempts failed. Agreement is measured against v2, never assumed; invalid/nonoptimal solutions are excluded.'}
    record('results/native_micom/settings.json',settings)
    controls=['baseline','no_external_lactate','block_D','block_L','block_both','block_acetate','B_no_lactate_uptake','no_lactate_block_D','no_lactate_block_L','no_lactate_block_both','no_lactate_block_acetate','no_lactate_B_no_uptake']
    rows=[];flux=[]
    for a in settings['abundances']:
        c=Community(pd.DataFrame({'id':['M','B'],'abundance':[a,1-a],'file':[str(x) for x in FILES.values()]}),solver='hybrid',max_exchange=100,progress=False)
        c.solver.problem.settings.update(optimality_tolerance=1e-7,primal_inf_tolerance=1e-8,dual_inf_tolerance=1e-8,iteration_limit=20000,time_limit=5,presolve=True,threads=1)
        med=pd.read_csv(BASE/'media/MJR8396C_gut_like_medium.tsv',sep='\t').set_index('compounds').maxFlux
        c.medium={ex.id:float(med[code]) for ex in c.exchanges for code in re.findall('cpd[0-9]+',ex.id) if code in med.index and med[code]>0}
        for r in c.reactions:
            if 'EX_cpd00007_' in r.id:r.bounds=(0,0)
        record(f'results/native_micom/medium_a{a}.json',c.medium)
        for ctrl in controls:
            with c:
                if ctrl.startswith('no_external') or ctrl.startswith('no_lactate_'):
                    for code in ['cpd00159','cpd00221']:
                        rid='EX_'+code+'_m'
                        if rid in c.reactions:c.reactions.get_by_id(rid).lower_bound=0
                base_ctrl=ctrl.replace('no_lactate_','')
                block={'block_D':['cpd00221'],'block_L':['cpd00159'],'block_both':['cpd00159','cpd00221'],'block_acetate':['cpd00029']}.get(base_ctrl,[])
                for code in block:
                    rid='EX_'+code+'_e0__B'
                    if rid in c.reactions:c.reactions.get_by_id(rid).upper_bound=0
                if ctrl in ['B_no_lactate_uptake','no_lactate_B_no_uptake']:
                    for code in ['cpd00159','cpd00221']:
                        rid='EX_'+code+'_e0__B'
                        if rid in c.reactions:c.reactions.get_by_id(rid).lower_bound=0
                gmax=c.slim_optimize()
                for f in settings['fractions']:
                    key={'M_abundance':a,'tradeoff':f,'control':ctrl}
                    try:
                        sol=c.cooperative_tradeoff(fraction=f,fluxes=True,pfba=True,atol=settings['atol'],rtol=settings['rtol'])
                        if sol is None:raise RuntimeError('Native MICOM returned None')
                        if sol.status!='optimal':raise RuntimeError('Nonoptimal native MICOM solution: '+sol.status)
                        growth={x:float(sol.members.loc[x,'growth_rate']) for x in ['M','B']}
                        row=key|{'status':sol.status,'community_max':gmax,'growth_M':growth['M'],'growth_B':growth['B'],'community_growth':sol.growth_rate,'weighted_growth':a*growth['M']+(1-a)*growth['B'],'regularizer':growth['M']**2+growth['B']**2}
                        sol.members.to_csv(out/f'members_a{a}_{ctrl}_f{f}.csv')
                        sol.fluxes.to_csv(out/f'fluxes_a{a}_{ctrl}_f{f}.csv')
                        for code in TARGETS:
                            rid='EX_'+code+'_e0';mid='EX_'+code+'_m'
                            M=float(sol.fluxes.loc['M'].get(rid,np.nan));B=float(sol.fluxes.loc['B'].get(rid,np.nan));pool=float(sol.fluxes.loc['medium'].get(mid,np.nan))
                            flux.append(key|{'compound':code,'M_specific':M,'B_specific':B,'M_weighted':a*M,'B_weighted':(1-a)*B,'community':pool,'accounting_residual':a*M+(1-a)*B-pool})
                        row['max_flux_accounting_residual']=max(abs(x['accounting_residual']) for x in flux[-len(TARGETS):] if np.isfinite(x['accounting_residual']))
                        if row['max_flux_accounting_residual']>1e-6 or row['community_growth']<f*gmax-1e-6:raise RuntimeError('Native solution failed mass/target growth validation')
                    except Exception as e:
                        row=key|{'status':'failed','error':repr(e)}
                        (ROOT/f'logs/micom_a{a}_{ctrl}_f{f}.txt').write_text(traceback.format_exc())
                    rows.append(row)
                    pd.DataFrame(rows).to_csv(out/'native_tradeoff_controls.csv',index=False)
                    pd.DataFrame(flux).to_csv(out/'native_taxon_community_fluxes.csv',index=False)
                print(a,ctrl,rows[-1],flush=True)
    old=pd.read_csv(BASE/'extended/results/community/tradeoff_abundance_controls.csv');new=pd.DataFrame(rows)
    merged=new.merge(old,on=['M_abundance','tradeoff','control'],suffixes=('_native','_v2'))
    for col in ['growth_M','growth_B','community_max']:merged[col+'_difference']=merged[col+'_native']-merged[col+'_v2']
    merged.to_csv(out/'native_vs_independent_v2.csv',index=False)
if __name__=='__main__':main()
