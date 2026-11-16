"""Complete native MICOM grid using explicit, certificate-audited QP adapters."""
from native_analysis import *
from micom import Community
from optlang.hybrid_interface import HybridProblem
from conditional_gradient_qp import solve_conditional_gradient
import traceback,time

OUT=ROOT/'results/native_micom_v4'
CONTROLS=['baseline','no_external_lactate','block_D','block_L','block_both','block_acetate','B_no_lactate_uptake','no_lactate_block_D','no_lactate_block_L','no_lactate_block_both','no_lactate_block_acetate','no_lactate_B_no_uptake']

def backend(self):
    solve_conditional_gradient(self)

def community(a):
    c=Community(pd.DataFrame({'id':['M','B'],'abundance':[a,1-a],'file':[str(x) for x in FILES.values()]}),solver='hybrid',max_exchange=100,progress=False)
    c.solver.problem.settings.update(primal_inf_tolerance=1e-8,dual_inf_tolerance=1e-8,optimality_tolerance=1e-9,time_limit=120,presolve=True,threads=1,verbose=0)
    med=pd.read_csv(BASE/'media/MJR8396C_gut_like_medium.tsv',sep='\t').set_index('compounds').maxFlux
    c.medium={ex.id:float(med[code]) for ex in c.exchanges for code in re.findall('cpd[0-9]+',ex.id) if code in med.index and med[code]>0}
    for r in c.reactions:
        if 'EX_cpd00007_' in r.id:r.bounds=(0,0)
    return c

def control(c,ctrl):
    if ctrl.startswith('no_external') or ctrl.startswith('no_lactate_'):
        for code in ['cpd00159','cpd00221']:c.reactions.get_by_id('EX_'+code+'_m').lower_bound=0
    base_ctrl=ctrl.replace('no_lactate_','')
    for code in {'block_D':['cpd00221'],'block_L':['cpd00159'],'block_both':['cpd00159','cpd00221'],'block_acetate':['cpd00029']}.get(base_ctrl,[]):
        c.reactions.get_by_id('EX_'+code+'_e0__B').upper_bound=0
    if ctrl in ['B_no_lactate_uptake','no_lactate_B_no_uptake']:
        for code in ['cpd00159','cpd00221']:c.reactions.get_by_id('EX_'+code+'_e0__B').lower_bound=0

def main():
    OUT.mkdir(exist_ok=True);HybridProblem.solve_osqp=backend
    rows=[];flux=[];certs=[]
    record('results/native_micom_v4/settings.json',{'engine':'native MICOM 0.39.1','MICOM_source_modified':False,'optlang_source_modified':False,'adapter':'explicit local certified fully corrective conditional-gradient QP backend; native full optlang matrix; HiGHS LP oracles and SciPy vertex QP; not a stock MICOM solver interface','native_cooperative_tradeoff':True,'native_pfba':True,'abundances':[.1,.25,.5,.75,.9],'tradeoffs':[.5,.8,.95,1.],'atol':1e-7,'rtol':1e-7,'solver_time_limit':120,'full_QP_certificate_and_full_network_checks':True,'prior_diagnostics':'native_micom_v4_multisolver_diagnostic; OSQP/Clarabel/PIQP/HiGHS direct-QP failures retained; no failed result accepted'})
    for a in [.1,.25,.5,.75,.9]:
        c=community(a)
        S=cobra.util.array.create_stoichiometric_matrix(c,array_type='lil').tocsr()
        for ctrl in CONTROLS:
            with c:
                control(c,ctrl);gmax=c.slim_optimize()
                for f in [.5,.8,.95,1.]:
                    key={'M_abundance':a,'tradeoff':f,'control':ctrl};start=time.monotonic()
                    try:
                        before=len(getattr(c.solver.problem,'qp_certificates',[]))
                        sol=c.cooperative_tradeoff(fraction=f,fluxes=True,pfba=True,atol=1e-7,rtol=1e-7)
                        attempts=c.solver.problem.qp_certificates[before:]
                        for k,cert in enumerate(attempts):
                            qvars=cert.pop('variables');trace=cert.pop('trace',None)
                            if trace:record(f'results/native_micom_v4/trace_a{a}_{ctrl}_f{f}_attempt{k}.json',trace)
                            certs.append(key|{'attempt':k,**cert})
                            pd.Series(qvars).to_csv(OUT/f'qp_variables_a{a}_{ctrl}_f{f}_attempt{k}.csv',header=['value'])
                        cert=attempts[-1]
                        if sol is None or sol.status!='optimal':raise RuntimeError('nonoptimal native MICOM result')
                        for col in ['max_primal_violation','max_stationarity_residual','max_complementarity_residual']:
                            if not np.isfinite(cert[col]) or cert[col]>1e-6:raise RuntimeError('QP certificate rejected: '+col+'='+str(cert[col]))
                        values=np.array([sol.fluxes.loc[r.community_id,r.global_id] for r in c.reactions])
                        residual=float(np.abs(S@values).max())
                        bounds=max(max(r.lower_bound-v,v-r.upper_bound,0) for r,v in zip(c.reactions,values))
                        gm,gb=map(float,sol.members.loc[['M','B'],'growth_rate'])
                        row=key|{'status':'optimal','community_max':gmax,'growth_M':gm,'growth_B':gb,'community_growth':sol.growth_rate,'weighted_growth':a*gm+(1-a)*gb,'regularizer':gm*gm+gb*gb,'mass_residual':residual,'bound_violation':bounds,'QP_backend':cert['backend'],'QP_objective':cert['objective'],'QP_global_gap':cert['global_suboptimality_bound'],'QP_feasibility':cert['max_primal_violation'],'QP_stationarity':cert['max_stationarity_residual'],'QP_complementarity':cert['max_complementarity_residual'],'elapsed_seconds':time.monotonic()-start}
                        if residual>1e-6 or bounds>1e-6 or sol.growth_rate<f*gmax-1e-6:raise RuntimeError('full pFBA feasibility rejected')
                        sol.fluxes.to_csv(OUT/f'fluxes_a{a}_{ctrl}_f{f}.csv');sol.members.to_csv(OUT/f'members_a{a}_{ctrl}_f{f}.csv')
                        for code in TARGETS:
                            m=float(sol.fluxes.loc['M'].get('EX_'+code+'_e0',np.nan));b=float(sol.fluxes.loc['B'].get('EX_'+code+'_e0',np.nan));pool=float(sol.fluxes.loc['medium'].get('EX_'+code+'_m',np.nan))
                            flux.append(key|{'compound':code,'M_specific':m,'B_specific':b,'M_weighted':a*m,'B_weighted':(1-a)*b,'community':pool,'accounting_residual':a*m+(1-a)*b-pool})
                    except Exception as e:
                        row=key|{'status':'failed','error':repr(e)}
                        (ROOT/f'logs/v4_micom_a{a}_{ctrl}_f{f}.txt').write_text(traceback.format_exc())
                    rows.append(row)
                    pd.DataFrame(rows).to_csv(OUT/'native_tradeoff_controls.csv',index=False)
                    pd.DataFrame(flux).to_csv(OUT/'native_taxon_community_fluxes.csv',index=False)
                    pd.DataFrame(certs).to_csv(OUT/'QP_certificates.csv',index=False)
                    print(key,row['status'],row.get('error',''),flush=True)
    old=pd.read_csv(BASE/'extended/results/community/tradeoff_abundance_controls.csv');new=pd.DataFrame(rows)
    merged=new.merge(old,on=['M_abundance','tradeoff','control'],suffixes=('_native','_v2'))
    for col in ['growth_M','growth_B','community_max']:merged[col+'_difference']=merged[col+'_native']-merged[col+'_v2']
    merged.to_csv(OUT/'native_vs_independent_v2.csv',index=False)

if __name__=='__main__':main()
