"""Diagnose comparator-medium and biomass limitations without gapfilling."""
from comparators import *

def main():
    audit=json.loads((ROOT/'results/comparators/structural_native_settings_audit.json').read_text())
    rows=[];masses=[]
    for source in audit:
        if 'AGORA' not in source['model']:continue
        m=load(path=ROOT/source['file']);mp={}
        for code,name in AGORA.items():
            for rid in ['EX_'+name+'_e_','EX_'+name+'(e)','EX_'+name+'[e]']:
                if rid in m.reactions:mp[code]=rid;break
        med=pd.read_csv(BASE/'media/MJR8396C_gut_like_medium.tsv',sep='\t')
        mapped={mp[r.compounds]:r.maxFlux for r in med.itertuples() if r.compounds in mp and r.maxFlux>0}
        m.medium=mapped;growth=m.slim_optimize();obj=list(cobra.util.solver.linear_reaction_coefficients(m))[0]
        b=m.reactions.get_by_id(source['biomass_audit']['synthesis_reaction']);atoms=collections.defaultdict(float);excluded=[]
        for met,value in b.metabolites.items():
            # Explicit bookkeeping placeholders produced from nothing are not
            # chemical polymer mass. Residual polymer X/R remains unresolved.
            if met.id in ['dnarep[c]','proteinsynth[c]','rnatrans[c]','biomass[c]']:
                excluded.append(met.id);continue
            for e,n in met.elements.items():atoms[e]-=value*n
        known=sum(MW.get(k,0)*v for k,v in atoms.items())/1000
        unknown={k:v for k,v in atoms.items() if k not in MW and abs(v)>1e-8}
        masses.append({'model':source['model'],'biomass_reaction':b.id,'known_element_net_mass_g_per_model_unit':known,'unresolved_elements':json.dumps(unknown),'bookkeeping_pseudometabolites_excluded':';'.join(excluded),'accepted_absolute_biomass_calibration':False,'reason':'Model formula/protonation/polymer conventions and experimental dry-weight calibration unresolved; mass accounting alone is not a growth calibration'})
        if growth>.1:
            rows.append({'model':source['model'],'base_growth':growth,'diagnostic_status':'already_exceeds_growth_target_0.1','extra_imports_total':0});continue
        obj.lower_bound=.1;penalties={}
        for r in m.exchanges:
            if r.id==mp.get('cpd00007'):continue
            if r.id in mapped:continue
            met,coef=next(iter(r.metabolites.items()))
            if coef<0:r.lower_bound=-1000;penalties[r.reverse_variable]=1.
            else:r.upper_bound=1000;penalties[r.forward_variable]=1.
        m.objective=m.problem.Objective(0,direction='min');m.objective.set_linear_coefficients(penalties)
        s=m.optimize()
        if s.status=='optimal':
            for r in m.exchanges:
                met,coef=next(iter(r.metabolites.items()));uptake=s.fluxes[r.id]*coef
                if r.id not in mapped and uptake>1e-7:
                    rows.append({'model':source['model'],'base_growth':growth,'diagnostic_status':s.status,'extra_imports_total':s.objective_value,'exchange':r.id,'metabolite':met.id,'name':met.name,'import_at_growth_0.1':uptake,'ModelSEED_annotation':met.annotation.get('seed.compound','')})
        else:rows.append({'model':source['model'],'base_growth':growth,'diagnostic_status':s.status})
        pd.DataFrame(rows).to_csv(ROOT/'results/comparators/v4_minimal_extra_imports_diagnostic.csv',index=False)
        print(source['model'],growth,s.status,flush=True)
    pd.DataFrame(rows).to_csv(ROOT/'results/comparators/v4_minimal_extra_imports_diagnostic.csv',index=False)
    pd.DataFrame(masses).to_csv(ROOT/'results/comparators/v4_biomass_formula_audit.csv',index=False)
    ime=load(path=ROOT/'inputs/comparators/iME375/Supplementary data 4.xml')
    record('results/comparators/v4_iME375_supplement_discrepancies.json',{'actual_biomass_synthesis':'R615','biomass_export':'R616 is reversible and written external-to-internal','growth_drain':'R617','deposited_biomass_weights':{x.id:v for x,v in ime.reactions.R615.metabolites.items()},'protein_synthesis_R610':ime.reactions.R610.reaction,'S1_serine_mmol_per_g_protein':.357,'deposited_R610_serine':-ime.reactions.R610.get_coefficient('LSER[c]'),'S1_leucine_mmol_per_g_protein':.492,'deposited_R610_leucine':-ime.reactions.R610.get_coefficient('LLEU[c]'),'decision':'Preserve deposited comparator; report source discrepancies, no silent repair or absolute-growth ranking'})
    record('results/comparators/v4_limits.json',{'diagnostic':'Minimum summed additional import at biomass flux 0.1; all unprovided exchanges except oxygen may open. This is a model-feasibility diagnostic, not an experimental auxotrophy or an adopted physiological medium.','harmonized_coculture':'Not promoted: the extra imports include amino acids, cofactors, dipeptides, fatty acids and sometimes other carbon sources; model biomass conventions also differ. Adding them selectively to force a comparative growth ranking would confound the experiment.','conclusion':'Comparator raw native/mapped simulations remain useful diagnostics; they do not establish absolute strain-specific growth or community equivalence.'})

if __name__=='__main__':main()
