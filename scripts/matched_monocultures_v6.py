"""Matched new strain GEMs: media, carbon, AA, maintenance, designs, FVA.

All predictions are computational comparisons under imposed bounds. Original
published NM10/BL7 models were unavailable; these newly generated GEMs are named
and pinned accordingly. Missing exchanges/targets are recorded, not invented.
"""
from native_analysis import *
import itertools,ast,argparse,gzip,faulthandler
from empirical_biomass_evidence import AA

OUT=ROOT/'results/matched_monocultures_v6'

def minimal_off(node):
    if isinstance(node,ast.Expression):return minimal_off(node.body)
    if isinstance(node,ast.Name):return [frozenset([node.id])]
    if isinstance(node,ast.BoolOp):
        parts=[minimal_off(x) for x in node.values]
        if isinstance(node.op,ast.And):return prune([s for p in parts for s in p])
        return prune([frozenset().union(*c) for c in itertools.product(*parts)])
    raise ValueError(ast.dump(node))

def prune(values):
    out=[]
    for s in sorted(set(values),key=lambda s:(len(s),sorted(s))):
        if not any(a<=s for a in out):out.append(s)
    return out

def main(name):
    faulthandler.dump_traceback_later(60,repeat=True)
    out=OUT/name;out.mkdir(parents=True,exist_ok=True)
    path=ROOT/'results/matched_reconstruction_v6'/name/(name+'.xml')
    m=load(path=path);bio=next(r for r in m.reactions if r.id=='bio1')
    rows=[];fvas=[];sens=[];aap=[];energy=[];prec=[];designs=[];envelopes=[]
    def collect(label,**more):
        row={'model':name,'scenario':label,**more,**phenotype(m)};rows.append(row);return row
    def fva(label,ids,method=None):
        for fraction in [.9,.95,1.]:
            try:
                f=flux_variability_analysis(m,reaction_list=ids,fraction_of_optimum=fraction,loopless=method,processes=1).reset_index(names='reaction')
                f['model']=name;f['scenario']=label;f['fraction']=fraction;f['loopless_method']=method or 'none';fvas.append(f)
            except Exception as e:fvas.append(pd.DataFrame([{'model':name,'scenario':label,'fraction':fraction,'loopless_method':method,'error':repr(e)}]))
    ids=['EX_'+c+'_e0' for c in TARGETS if 'EX_'+c+'_e0' in m.reactions]
    for med in ['gut_like','carbon_rich','AA_rich']:
        missing=medium(m,med);collect(med,unmapped_medium=';'.join(missing));fva(med,ids)
    print(name,'media FVA complete',flush=True)
    # Carbon source controls retain the same specified six-AA/vitamin background.
    medium(m)
    for c in ['cpd00027','cpd00159','cpd00221','cpd00029']:
        if 'EX_'+c+'_e0' in m.reactions:m.reactions.get_by_id('EX_'+c+'_e0').lower_bound=0
    bg=collect('carbon_background')['growth']
    panel=pd.read_csv(BASE/'tables/Table_S11_carbon_normalized_substrate_screen.csv')
    for p in panel.itertuples():
        with m:
            rid='EX_'+p.compound+'_e0'
            if rid not in m.reactions:
                rows.append({'model':name,'scenario':'carbon60_'+p.substrate,'compound':p.compound,'status':'unmapped_exchange','background_growth':bg});continue
            m.reactions.get_by_id(rid).lower_bound=-60/float(p.carbon_atoms)
            rr=collect('carbon60_'+p.substrate,compound=p.compound,carbon_supply=60,background_growth=bg)
            rr['incremental_growth']=rr['growth']-bg
            if p.compound in ['cpd00027','cpd00159','cpd00221','cpd00029']:fva('carbon60_'+p.substrate,ids)
    print(name,'carbon screen complete',flush=True)
    # Acetate dependence at fixed glucose and no supplied lactate.
    medium(m)
    for c in ['cpd00159','cpd00221']:
        if 'EX_'+c+'_e0' in m.reactions:m.reactions.get_by_id('EX_'+c+'_e0').lower_bound=0
    for cap in [0,5,10,20,25]:
        with m:
            if 'EX_cpd00029_e0' in m.reactions:m.reactions.EX_cpd00029_e0.lower_bound=-cap
            collect('acetate_cosubstrate',acetate_capacity=cap)
    # All twenty individual AA additions and all 190 pairs, including canonical
    # Stickland donor/acceptor substrates; flux outputs accompany growth.
    medium(m)
    for c in AA.values():
        rid='EX_'+c+'_e0'
        if rid in m.reactions:m.reactions.get_by_id(rid).lower_bound=0
    base=phenotype(m);single={}
    for aa,c in AA.items():
        rid='EX_'+c+'_e0'
        with m:
            if rid not in m.reactions:
                aap.append({'model':name,'AA1':aa,'status':'unmapped_exchange'});continue
            m.reactions.get_by_id(rid).lower_bound=-1
            row={'model':name,'AA1':aa,'AA2':'','background_growth':base['growth'],'uptake_capacity_each':1,**phenotype(m)}
            single[aa]=row['growth'];aap.append(row)
    for a,b in itertools.combinations(AA,2):
        with m:
            rids=['EX_'+AA[x]+'_e0' for x in [a,b]]
            if not all(r in m.reactions for r in rids):
                aap.append({'model':name,'AA1':a,'AA2':b,'status':'unmapped_exchange'});continue
            for r in rids:m.reactions.get_by_id(r).lower_bound=-1
            row={'model':name,'AA1':a,'AA2':b,'background_growth':base['growth'],'uptake_capacity_each':1,**phenotype(m)}
            row['interaction_effect']=row['growth']-single[a]-single[b]+base['growth'];aap.append(row)
    print(name,'20 AA singles and 190 pairs complete',flush=True)
    medium(m);gam=-bio.get_coefficient('cpd00008_c0') # ADP is a product; use its positive coefficient.
    gam=bio.get_coefficient('cpd00008_c0')
    for med,scale,GAM,NGAM in itertools.product(['gut_like','carbon_rich','AA_rich'],[.5,1,1.5],[20,40,60],[0,1,2,5]):
        with m:
            medium(m,med)
            for ex in m.exchanges:
                met=next(iter(ex.metabolites))
                if ex.lower_bound<0 and met.elements.get('C',0)>0 and met.id!='cpd00011_e0':ex.lower_bound*=scale
            change=GAM-gam
            bio.add_metabolites({m.metabolites.get_by_id(c+'_c0'):v*change for c,v in {'cpd00002':-1,'cpd00001':-1,'cpd00008':1,'cpd00009':1,'cpd00067':1}.items()})
            m.reactions.rxn00062_c0.lower_bound=NGAM
            sens.append({'model':name,'medium':med,'organic_uptake_scale':scale,'GAM_assumed':GAM,'NGAM_assumed':NGAM,**phenotype(m)})
    print(name,'maintenance grid complete',flush=True)
    medium(m)
    for met,v in list(bio.metabolites.items()):
        if v>=0:continue
        with m:
            dm=m.add_boundary(met,type='demand',reaction_id='V6_PRECURSOR_DEMAND');m.objective=dm
            prec.append({'model':name,'metabolite':met.id,'max_demand':m.slim_optimize(error_value=np.nan),'status':m.solver.status})
    for condition in ['closed_boundaries','inorganics_only']:
        with m:
            for r in m.boundary:r.bounds=(0,0)
            if condition=='inorganics_only':
                for c in ['cpd00001','cpd00009','cpd00011','cpd00013','cpd00067','cpd00048']:
                    rid='EX_'+c+'_e0'
                    if rid in m.reactions:m.reactions.get_by_id(rid).bounds=(-100,100)
            atp=m.reactions.rxn00062_c0;atp.bounds=(0,1000);m.objective=atp
            energy.append({'model':name,'condition':condition,'ATP_max':m.slim_optimize(error_value=np.nan),'status':m.solver.status})
    # Same chemical intervention targets, with native fresh GPRs and global
    # minimum gene-set implementations where every target has an assigned GPR.
    history=pd.read_csv(BASE/'extended/results/strain_design/historical_reaction_deletion_retests.csv').fillna('')
    medium(m)
    for row in history.itertuples():
        rids=row.target_reactions.split(';') if row.target_reactions else []
        if any(r not in m.reactions for r in rids):
            designs.append({'model':name,'design':row.design,'level':'reaction','status':'unrepresented_target','missing':';'.join(r for r in rids if r not in m.reactions)});continue
        with m:
            for rid in rids:m.reactions.get_by_id(rid).knock_out()
            ph=phenotype(m);designs.append({'model':name,'design':row.design,'level':'reaction','targets':';'.join(rids),**ph})
            if 'EX_cpd00211_e0' in m.reactions and ph.get('growth',0)>1e-8:
                for frac in [0,.25,.5,.75,.9,.95,1.]:
                    with m:
                        obj=next(iter(cobra.util.solver.linear_reaction_coefficients(m)))
                        obj.lower_bound=max(obj.lower_bound,frac*ph['growth']);m.objective=m.reactions.EX_cpd00211_e0
                        m.objective_direction='min';lo=m.slim_optimize(error_value=np.nan);m.objective_direction='max';hi=m.slim_optimize(error_value=np.nan)
                        envelopes.append({'model':name,'design':row.design,'growth_fraction':frac,'growth_floor':frac*ph['growth'],'butyrate_min':lo,'butyrate_max':hi})
        if any(not m.reactions.get_by_id(rid).gene_reaction_rule for rid in rids):
            designs.append({'model':name,'design':row.design,'level':'gene','status':'missing_GPR_cannot_implement'});continue
        choices=[minimal_off(m.reactions.get_by_id(rid).gpr.body) for rid in rids]
        sets=prune([frozenset().union(*cs) for cs in itertools.product(*choices)]) if choices else [frozenset()]
        for gs in sets:
            with m:
                cobra.manipulation.knock_out_model_genes(m,list(gs))
                designs.append({'model':name,'design':row.design,'level':'gene','genes':';'.join(sorted(gs)),'targets':';'.join(rids),**phenotype(m)})
    # Loopless calculations run separately with a bounded native fastSNP
    # solver. Keep the complete ordinary results even if that optional
    # numerical step stalls on a newly reconstructed network.
    for file,data in [('phenotypes',rows),('maintenance_uptake_grid',sens),('AA_singles_and_pairs',aap),('ATP_cycle_checks',energy),('precursor_producibility',prec),('intervention_implementations',designs),('growth_product_envelopes',envelopes)]:pd.DataFrame(data).to_csv(out/(file+'.csv'),index=False)
    pd.concat(fvas).to_csv(out/'exchange_FVA.csv',index=False)
    draft=load(path=path.with_name(name+'-draft.xml'));added=sorted(set(m.reactions.list_attr('id'))-set(draft.reactions.list_attr('id')))
    pd.DataFrame([{'reaction':rid,'name':m.reactions.get_by_id(rid).name,'GPR':m.reactions.get_by_id(rid).gene_reaction_rule} for rid in added]).to_csv(out/'gapfill_additions.csv',index=False)
    record(str((out/'model_audit.json').relative_to(ROOT)),{'model':name,'new_automated_reconstruction':True,'sha256':sha(path),'reactions':len(m.reactions),'metabolites':len(m.metabolites),'genes':len(m.genes),'draft_reactions':len(draft.reactions),'added_reactions':len(added),'biomass_mass_g_per_model_unit':mass(bio.metabolites),'assumed_GAM':gam,'no_measured_target_biomass':True,'AA_pairs':190,'AA_single_panel':20,'maintenance_grid':108})
    faulthandler.cancel_dump_traceback_later()
    print(name,'matched monoculture complete',len(rows),len(sens),len(aap),len(designs),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--strain',required=True,choices=['MJR8396C','NM10','BL7','B22L']);a=p.parse_args();main(a.strain)
