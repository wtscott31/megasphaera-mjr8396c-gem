"""Native MICOM LP maxima and feasibility/pFBA at supplied v2 growth rates.

This validates MICOM's native community construction independently of its
numerically unstable full-network OSQP tradeoff stage. Fixed rates are explicitly
imported from v2; these rows are not native cooperative-QP optimizations.
"""
from native_micom import *
def main():
    old=pd.read_csv(BASE/'extended/results/community/tradeoff_abundance_controls.csv')
    out=ROOT/'results/native_micom_linear';out.mkdir(exist_ok=True)
    rows=[];fixed=[];flux=[];fv=[]
    controls=old.control.unique()
    for a in [.1,.25,.5,.75,.9]:
        c=Community(pd.DataFrame({'id':['M','B'],'abundance':[a,1-a],'file':[str(x) for x in FILES.values()]}),solver='glpk',max_exchange=100,progress=False)
        c.solver.configuration.tolerances.feasibility=1e-8
        med=pd.read_csv(BASE/'media/MJR8396C_gut_like_medium.tsv',sep='\t').set_index('compounds').maxFlux
        c.medium={ex.id:float(med[code]) for ex in c.exchanges for code in re.findall('cpd[0-9]+',ex.id) if code in med.index and med[code]>0}
        for r in c.reactions:
            if 'EX_cpd00007_' in r.id:r.bounds=(0,0)
        for ctrl in controls:
            with c:
                if ctrl.startswith('no_external') or ctrl.startswith('no_lactate_'):
                    for code in ['cpd00159','cpd00221']:
                        rid='EX_'+code+'_m'
                        if rid in c.reactions:c.reactions.get_by_id(rid).lower_bound=0
                base_ctrl=ctrl.replace('no_lactate_','')
                for code in {'block_D':['cpd00221'],'block_L':['cpd00159'],'block_both':['cpd00159','cpd00221'],'block_acetate':['cpd00029']}.get(base_ctrl,[]):
                    rid='EX_'+code+'_e0__B'
                    if rid in c.reactions:c.reactions.get_by_id(rid).upper_bound=0
                if ctrl in ['B_no_lactate_uptake','no_lactate_B_no_uptake']:
                    for code in ['cpd00159','cpd00221']:
                        rid='EX_'+code+'_e0__B'
                        if rid in c.reactions:c.reactions.get_by_id(rid).lower_bound=0
                maximum=c.slim_optimize();v2=old[(old.M_abundance==a)&(old.control==ctrl)].iloc[0]
                rows.append({'M_abundance':a,'control':ctrl,'native_LP_community_max':maximum,'v2_community_max':v2.community_max,'difference':maximum-v2.community_max,'status':c.solver.status})
                if a!=.5:continue
                for point in old[(old.M_abundance==a)&(old.control==ctrl)].itertuples():
                    with c:
                        for taxon,g in [('M',point.growth_M),('B',point.growth_B)]:
                            const=c.constraints['objective_'+taxon];const.ub=None;const.lb=max(0,g-1e-8);const.ub=g+1e-8
                        s=c.optimize(fluxes=True,pfba=True,raise_error=True,atol=1e-7,rtol=1e-7)
                        if s is None or s.status!='optimal':raise RuntimeError('Native fixed-growth LP failed')
                        key={'M_abundance':a,'control':ctrl,'tradeoff':point.tradeoff,'method':'native_MICOM_LP_pFBA_at_v2_supplied_member_growth'}
                        fixed.append(key|{'status':s.status,'growth_M':float(s.members.loc['M','growth_rate']),'growth_B':float(s.members.loc['B','growth_rate']),'community_growth':s.growth_rate})
                        s.fluxes.to_csv(out/f'fluxes_{ctrl}_f{point.tradeoff}.csv')
                        for code in TARGETS:
                            rid='EX_'+code+'_e0';mid='EX_'+code+'_m'
                            M=float(s.fluxes.loc['M'].get(rid,np.nan));B=float(s.fluxes.loc['B'].get(rid,np.nan));pool=float(s.fluxes.loc['medium'].get(mid,np.nan))
                            flux.append(key|{'compound':code,'M_specific':M,'B_specific':B,'M_weighted':a*M,'B_weighted':(1-a)*B,'community':pool,'accounting_residual':a*M+(1-a)*B-pool})
                        if point.tradeoff==.95 and ctrl in ['baseline','no_external_lactate','block_both','no_lactate_block_both']:
                            ids=[f'EX_{code}_e0__{taxon}' for code in ['cpd00159','cpd00221','cpd00211','cpd00029','cpd00141'] for taxon in ['M','B'] if f'EX_{code}_e0__{taxon}' in c.reactions]
                            f=flux_variability_analysis(c,reaction_list=ids,fraction_of_optimum=1.,processes=1).reset_index(names='reaction');f['control']=ctrl;fv.append(f)
                    # Raw optlang bound setters do not register a COBRA context
                    # undo. Restore explicitly before the next control's LP.
                    for taxon in ['M','B']:
                        const=c.constraints['objective_'+taxon];const.ub=None;const.lb=0
            pd.DataFrame(rows).to_csv(out/'native_LP_maxima_vs_v2.csv',index=False)
            pd.DataFrame(fixed).to_csv(out/'native_fixed_growth_validation.csv',index=False)
            pd.DataFrame(flux).to_csv(out/'native_fixed_growth_exchange_accounting.csv',index=False)
            if fv:pd.concat(fv).to_csv(out/'native_fixed_growth_exchange_FVA.csv',index=False)
        print('Native MICOM LP abundance',a,'complete',flush=True)
    pd.DataFrame(rows).to_csv(out/'native_LP_maxima_vs_v2.csv',index=False)
if __name__=='__main__':main()
