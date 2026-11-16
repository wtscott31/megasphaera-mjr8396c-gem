"""Curated deposited models: native and mapped anaerobic media, no gapfilling."""
from native_analysis import *
AGORA={'cpd00027':'glc_D','cpd00159':'lac_L','cpd00221':'lac_D','cpd00029':'ac','cpd00211':'but','cpd00141':'ppa','cpd00047':'for','cpd00036':'succ','cpd00013':'nh4','cpd00239':'h2s','cpd00324':'methanethiol','cpd00011':'co2','cpd11640':'h2','cpd00067':'h','cpd00001':'h2o','cpd00007':'o2','cpd00009':'pi','cpd00048':'so4','cpd00971':'na1','cpd00205':'k','cpd00254':'mg2','cpd00063':'ca2','cpd00099':'cl','cpd10515':'fe2','cpd10516':'fe3','cpd00030':'mn2','cpd00034':'zn2','cpd00058':'cu2','cpd00149':'cobalt2','cpd00305':'thm','cpd00220':'ribflv','cpd00133':'ncam','cpd00218':'nac','cpd00644':'pnto_R','cpd00104':'btn','cpd00041':'asp_L','cpd00033':'gly','cpd00039':'lys_L','cpd00060':'met_L','cpd00084':'cys_L','cpd00119':'his_L'}
IME={'cpd00027':'R001','cpd00221':'R018','cpd00029':'R042','cpd00211':'R043','cpd00141':'R044','cpd00047':'R017','cpd00013':'R014','cpd00011':'R010','cpd11640':'R011','cpd00001':'R012','cpd00009':'R021','cpd00048':'R020','cpd00033':'R022','cpd00060':'R029','cpd00084':'R030','cpd00041':'R031','cpd00119':'R036','cpd00039':'R040','cpd00218':'R049'}
def biomassaudit(m,obj):
    r=m.reactions.get_by_id(obj)
    # Trace drains and biomass-only transport to the actual synthesis reaction.
    visited=set()
    for _ in range(8):
        if r.id in visited:break
        visited.add(r.id)
        nonbio=[met for met in r.metabolites if 'biomass' not in (met.id+' '+met.name).lower()]
        if nonbio:break
        consumed=[met for met,v in r.metabolites.items() if v<0]
        if len(consumed)!=1:break
        met=consumed[0];candidates=[x for x in met.reactions if x.id not in visited and x.get_coefficient(met)>0]
        if len(candidates)!=1:break
        r=candidates[0]
    missing=[];atoms=collections.defaultdict(float)
    for met,coef in r.metabolites.items():
        if 'biomass' in (met.id+' '+met.name).lower():continue
        if not met.formula:missing.append(met.id);continue
        for k,v in met.elements.items():atoms[k]-=coef*v
    unresolved={k:v for k,v in atoms.items() if k not in MW and abs(v)>1e-8}
    value=sum(MW.get(k,0)*v for k,v in atoms.items())/1000 if not missing and not unresolved else None
    if value is not None and value<=1e-6:value=None
    return {'synthesis_reaction':r.id,'net_mass_model_formulas_g_per_biomass':value,'missing_formulas':missing,'unresolved_elements':unresolved,'interpretation':'Model formula convention; not empirical calibration. Unknown mass prevents absolute-growth equivalence claims.'}
def main():
    files={'MJR8396C':FILES['M'],'iME375':ROOT/'inputs/comparators/iME375/Supplementary data 4.xml'}
    files.update({'AGORA2_fixed_'+p.stem:p for p in (ROOT/'inputs/comparators/AGORA2_fixed').glob('*.xml')})
    for p in (ROOT/'inputs/comparators/encoding_repaired').glob('*.xml'):files['AGORA2_fixed_'+p.stem]=p
    out=ROOT/'results/comparators';out.mkdir(exist_ok=True);rows=[];maps=[];audit=[];fvas=[]
    for name,p in files.items():
        m=cobra.io.read_sbml_model(str(p));m.solver='glpk'
        original_obj=str(m.objective.expression);native_bounds={r.id:list(r.bounds) for r in m.exchanges}
        if name=='iME375':mapping={k:v for k,v in IME.items() if v in m.reactions};biomass='R617'
        elif name=='MJR8396C':mapping={k:'EX_'+k+'_e0' for k in AGORA if 'EX_'+k+'_e0' in m.reactions};biomass='EX_cpd11416_c0'
        else:
            mapping={}
            for k,v in AGORA.items():
                for rid in ['EX_'+v+'_e_','EX_'+v+'(e)','EX_'+v+'[e]']:
                    if rid in m.reactions:mapping[k]=rid;break
            biomass=list(cobra.util.solver.linear_reaction_coefficients(m))[0].id
        bd=biomassaudit(m,biomass)
        maintenance=[{'id':r.id,'bounds':r.bounds,'reaction':r.reaction} for r in m.reactions if r.id=='ATPM' or r.id=='rxn00062_c0' or 'atp maintenance' in r.name.lower()]
        audit.append({'model':name,'file':str(p.relative_to(ROOT)),'sha256':sha(p),'reactions':len(m.reactions),'metabolites':len(m.metabolites),'compartments':m.compartments,'deposited_objective':original_obj,'growth_objective':biomass,'biomass_audit':bd,'maintenance':maintenance,'native_exchange_bounds':native_bounds,'notes':'iME375 deposited objective is hexanoate, not growth; native objective retained as a separate row. AGORA deposited all-open exchanges are not a physiological published diet.'})
        def result(label,objective=None):
            if objective:m.objective=objective
            s=m.optimize();row={'model':name,'scenario':label,'objective':str(m.objective.expression),'status':s.status,'objective_value':s.objective_value,'growth_flux':float(s.fluxes[biomass]) if s.status=='optimal' else np.nan,'biomass_mass_g_per_unit':bd['net_mass_model_formulas_g_per_biomass']}
            if s.status=='optimal':
                try:
                    pf=pfba(m,fraction_of_optimum=1-1e-8)
                    row.update({k:float(pf.fluxes.get(v,np.nan)) for k,v in mapping.items() if k in TARGETS})
                    if s.objective_value>1e-8 and label in ['harmonized_gut','carbon60_glucose','carbon60_L_lactate_acetate5','carbon60_D_lactate_acetate5']:
                        ids=[mapping[x] for x in ['cpd00211','cpd00141','cpd00029','cpd00159','cpd00221'] if x in mapping]
                        f=flux_variability_analysis(m,reaction_list=ids,fraction_of_optimum=.95,processes=1).reset_index(names='reaction');f['model']=name;f['scenario']=label;fvas.append(f)
                except Exception as e:row['pFBA_error']=repr(e)
            rows.append(row)
        result('native_as_deposited_objective')
        with m:
            for code in ['cpd00007']:
                if code in mapping:m.reactions.get_by_id(mapping[code]).lower_bound=0
            result('native_bounds_growth_objective_anaerobic',biomass)
        med=pd.read_csv(BASE/'media/MJR8396C_gut_like_medium.tsv',sep='\t')
        mapped={}
        for r in med.itertuples():
            rid=mapping.get(r.compounds)
            maps.append({'model':name,'ModelSEED':r.compounds,'name':r.name,'uptake_capacity_mmol_per_model_gDW_h':r.maxFlux,'mapped_exchange':rid,'status':'mapped' if rid else 'missing_no_substitute','exchange_sign':next(iter(m.reactions.get_by_id(rid).metabolites.values())) if rid else None,'metabolite_name':next(iter(m.reactions.get_by_id(rid).metabolites)).name if rid else None})
            if rid and r.maxFlux>0:mapped[rid]=float(r.maxFlux)
        with m:
            m.medium=mapped;m.objective=biomass
            result('harmonized_gut')
            for substrate,carbon in [('glucose',6),('L_lactate',3),('D_lactate',3)]:
                code={'glucose':'cpd00027','L_lactate':'cpd00159','D_lactate':'cpd00221'}[substrate]
                for acetate in [0,5]:
                    label=f'carbon60_{substrate}'+('_acetate5' if acetate else '')
                    with m:
                        for c in ['cpd00027','cpd00159','cpd00221','cpd00029']:
                            if c in mapping:m.reactions.get_by_id(mapping[c]).lower_bound=0
                        if code not in mapping:
                            rows.append({'model':name,'scenario':label,'status':'substrate_exchange_unmapped','growth_flux':np.nan});continue
                        # Allow secretion when a deposited exchange was zeroed,
                        # recorded as an explicit scenario bound, not a model edit.
                        ex=m.reactions.get_by_id(mapping[code]);ex.lower_bound=-60/carbon
                        if acetate and 'cpd00029' in mapping:m.reactions.get_by_id(mapping['cpd00029']).lower_bound=-acetate
                        result(label)
        pd.DataFrame(rows).to_csv(out/'curated_comparator_simulations.csv',index=False)
        pd.DataFrame(maps).to_csv(out/'medium_mapping_ledger.csv',index=False)
        record('results/comparators/structural_native_settings_audit.json',audit)
        if fvas:pd.concat(fvas).to_csv(out/'selected_exchange_FVA_95percent.csv',index=False)
        print(name,'completed',flush=True)
    record('results/comparators/method_limits.json',{'carbon_supply':'60 mmol main-substrate C per model gDW per hour; six template amino acids and mapped vitamins retained equally where represented; acetate sensitivity supplies an additional 10 mmol C. Not a single-substrate defined medium.','missing_components':'Every unmapped component is listed; no additional nutrient, transporter, gapfill, or reaction was invented. Zero/infeasible growth is retained.','absolute_growth':'Compare only if biomass mass normalization is supported; iME375 has missing formulas and a product-oriented deposited objective. Raw fluxes are not proof of equivalent gDW basis.','coculture':'Curated-comparator cocultures not promoted: incomplete common nutrient/transport mapping and biomass normalization require resolution first. Native two-member deposited-model controls are separate.','strain':'Lee et al. 2020 section 2.1 specifies DSM20460; original HE576794 genome GCA_000283495.2 distinguished from later same-strain GCA_003010495.1.'})
if __name__=='__main__':main()
