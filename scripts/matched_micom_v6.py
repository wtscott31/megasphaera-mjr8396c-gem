"""New matched native MICOM strain comparisons plus 22L curation controls."""
from native_micom_complete import *
import argparse

def build_pair(mfile,bfile,a):
    c=Community(pd.DataFrame({'id':['M','B'],'abundance':[a,1-a],'file':[str(mfile),str(bfile)]}),solver='hybrid',max_exchange=100,progress=False)
    c.solver.problem.settings.update(primal_inf_tolerance=1e-8,dual_inf_tolerance=1e-8,optimality_tolerance=1e-9,time_limit=120,presolve=True,threads=1,verbose=0)
    med=pd.read_csv(BASE/'media/MJR8396C_gut_like_medium.tsv',sep='\t').set_index('compounds').maxFlux
    c.medium={ex.id:float(med[code]) for ex in c.exchanges for code in re.findall('cpd[0-9]+',ex.id) if code in med.index and med[code]>0}
    for r in c.reactions:
        if 'EX_cpd00007_' in r.id:r.bounds=(0,0)
    return c

def safe_control(c,ctrl):
    missing=[]
    if ctrl.startswith('no_external') or ctrl.startswith('no_lactate_'):
        for code in ['cpd00159','cpd00221']:
            rid='EX_'+code+'_m'
            if rid in c.reactions:c.reactions.get_by_id(rid).lower_bound=0
    base_ctrl=ctrl.replace('no_lactate_','')
    for code in {'block_D':['cpd00221'],'block_L':['cpd00159'],'block_both':['cpd00159','cpd00221'],'block_acetate':['cpd00029']}.get(base_ctrl,[]):
        rid='EX_'+code+'_e0__B'
        if rid in c.reactions:c.reactions.get_by_id(rid).upper_bound=0
        else:missing.append(rid)
    if ctrl in ['B_no_lactate_uptake','no_lactate_B_no_uptake']:
        for code in ['cpd00159','cpd00221']:
            rid='EX_'+code+'_e0__B'
            if rid in c.reactions:c.reactions.get_by_id(rid).lower_bound=0
            else:missing.append(rid)
    return missing

def pairs():
    d=ROOT/'results/matched_reconstruction_v6'
    x={n:(d/n/(n+'.xml'),d/'B22L/B22L.xml') for n in ['MJR8396C','NM10','BL7']}
    x['GPRcurated_original22L']=(ROOT/'models/MJR8396C_sequence_supported_GPR_curated_v6.xml',FILES['B'])
    x['GPRcurated_filtered22L']=(ROOT/'models/MJR8396C_sequence_supported_GPR_curated_v6.xml',ROOT/'models/B22L_D_LDH_unsupported_gapfill_disabled_v6.xml')
    for n in ['MJR8396C','NM10','BL7']:
        x[n+'_matched_curated_sensitivity']=(ROOT/'models'/(n+'_matched_sequence_curated_v6.xml'),ROOT/'models/B22L_matched_sequence_curated_v6.xml')
    return x

def run_pair(name,mfile,bfile):
    out=ROOT/'results/matched_micom_v6'/name;out.mkdir(parents=True,exist_ok=True)
    assert mfile.exists() and bfile.exists()
    HybridProblem.solve_osqp=backend
    rows=[];flux=[];certs=[];fv=[];mono=[]
    record(str((out/'settings.json').relative_to(ROOT)),{'native_MICOM':'0.39.1','local_QP_adapter':'conditional_gradient_qp.py; full native matrix; no MICOM source edits',
          'M':str(mfile.relative_to(ROOT)),'M_sha256':sha(mfile),'B':str(bfile.relative_to(ROOT)),'B_sha256':sha(bfile),
          'medium_basis':'single shared pool mmol / community gDW / h; abundance-weighted taxon exchange accounting',
          'controls':CONTROLS,'abundances':[.1,.25,.5,.75,.9],'tradeoffs':[.5,.8,.95,1.],
          'reconstruction_scope':'MJR8396C/NM10/BL7 pairs use NEW uniform automated reconstructions; GPRcurated pairs derive from deposited models'})
    # Effective per-species maximum supply is pooled bound / abundance when
    # compared to the species alone; total-specific bounds are not duplicated.
    for a in [.1,.25,.5,.75,.9]:
        for label,p,ab in [('M',mfile,a),('B',bfile,1-a)]:
            m=load(path=p);medium(m)
            m.medium={rid:v/ab for rid,v in m.medium.items()}
            for ctrl in ['baseline','no_external_lactate']:
                with m:
                    if ctrl=='no_external_lactate':
                        for code in ['cpd00159','cpd00221']:
                            rid='EX_'+code+'_e0'
                            if rid in m.reactions:m.reactions.get_by_id(rid).lower_bound=0
                    mono.append({'pair':name,'member':label,'abundance':ab,'M_abundance':a,'control':ctrl,**phenotype(m)})
        c=build_pair(mfile,bfile,a);S=cobra.util.array.create_stoichiometric_matrix(c,array_type='lil').tocsr()
        for ctrl in CONTROLS:
            with c:
                missing=safe_control(c,ctrl);gmax=c.slim_optimize()
                for f in [.5,.8,.95,1.]:
                    key={'pair':name,'M_abundance':a,'tradeoff':f,'control':ctrl};start=time.monotonic()
                    try:
                        before=len(getattr(c.solver.problem,'qp_certificates',[]))
                        sol=c.cooperative_tradeoff(fraction=f,fluxes=True,pfba=True,atol=1e-7,rtol=1e-7)
                        attempts=c.solver.problem.qp_certificates[before:]
                        cert=attempts[-1]
                        if sol is None or sol.status!='optimal':raise RuntimeError('Nonoptimal MICOM solution')
                        for col in ['max_primal_violation','max_stationarity_residual','max_complementarity_residual']:
                            if not np.isfinite(cert[col]) or cert[col]>1e-6:raise RuntimeError('QP certificate rejected '+col)
                        values=np.array([sol.fluxes.loc[r.community_id,r.global_id] for r in c.reactions])
                        residual=float(abs(S@values).max());bounds=max(max(r.lower_bound-v,v-r.upper_bound,0) for r,v in zip(c.reactions,values))
                        gm,gb=map(float,sol.members.loc[['M','B'],'growth_rate'])
                        if residual>1e-6 or bounds>1e-6 or sol.growth_rate<f*gmax-1e-6:raise RuntimeError('Full network feasibility rejected')
                        row=key|{'status':'optimal','community_max':gmax,'growth_M':gm,'growth_B':gb,'community_growth':sol.growth_rate,'weighted_growth':a*gm+(1-a)*gb,
                                 'mass_residual':residual,'bound_violation':bounds,'QP_global_gap':cert['global_suboptimality_bound'],
                                 'QP_feasibility':cert['max_primal_violation'],'QP_stationarity':cert['max_stationarity_residual'],
                                 'QP_complementarity':cert['max_complementarity_residual'],'elapsed_seconds':time.monotonic()-start,'unrepresented_control_targets':';'.join(missing)}
                        for k,ce in enumerate(attempts):
                            qvars=ce.pop('variables');trace=ce.pop('trace',None)
                            certs.append(key|{'attempt':k,**ce})
                            pd.Series(qvars).to_csv(out/f'qp_variables_a{a}_{ctrl}_f{f}_attempt{k}.csv',header=['value'])
                        sol.fluxes.to_csv(out/f'fluxes_a{a}_{ctrl}_f{f}.csv')
                        for code in TARGETS:
                            vals=[sol.fluxes.loc[taxon].get(rid,0) for taxon,rid in [('M','EX_'+code+'_e0'),('B','EX_'+code+'_e0'),('medium','EX_'+code+'_m')]]
                            mv,bv,pool=[float(v) if pd.notna(v) else 0. for v in vals]
                            flux.append(key|{'compound':code,'M_specific':mv,'B_specific':bv,'M_weighted':a*mv,'B_weighted':(1-a)*bv,'community':pool,'accounting_residual':a*mv+(1-a)*bv-pool})
                    except Exception as e:
                        row=key|{'status':'failed','error':repr(e),'elapsed_seconds':time.monotonic()-start}
                        (out/f'error_a{a}_{ctrl}_f{f}.txt').write_text(traceback.format_exc())
                    rows.append(row)
                # Exchange FVA on the full community at >=95% maximal growth;
                # this measures feasible diversity, not uncertainty in the QP point.
                ids=[rid for code in ['cpd00159','cpd00221','cpd00211','cpd00141','cpd00029'] for rid in ['EX_'+code+'_e0__M','EX_'+code+'_e0__B'] if rid in c.reactions]
                try:
                    r=flux_variability_analysis(c,reaction_list=ids,fraction_of_optimum=.95,processes=1).reset_index(names='reaction')
                    r['pair']=name;r['M_abundance']=a;r['control']=ctrl;r['objective_fraction']=.95;fv.append(r)
                except Exception as e:
                    fv.append(pd.DataFrame([{'pair':name,'M_abundance':a,'control':ctrl,'error':repr(e)}]))
            pd.DataFrame(rows).to_csv(out/'native_tradeoff_controls.csv',index=False)
            pd.DataFrame(flux).to_csv(out/'native_taxon_community_fluxes.csv',index=False)
            pd.DataFrame(certs).to_csv(out/'QP_certificates.csv',index=False)
            pd.DataFrame(mono).to_csv(out/'effective_supply_monocultures.csv',index=False)
            if fv:pd.concat(fv).to_csv(out/'exchange_FVA95.csv',index=False)
        print(name,'abundance',a,'completed',sum(r['status']=='optimal' for r in rows),'/',len(rows),flush=True)
    return rows

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--pair',choices=list(pairs()),required=True);a=p.parse_args()
    run_pair(a.pair,*pairs()[a.pair])
