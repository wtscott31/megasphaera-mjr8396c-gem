"""Distinguish incorporated biomass moieties from recycled carrier pools."""
from native_analysis import *

def main():
    rows=[];ledger=[]
    for name in ['MJR8396C','NM10','BL7']:
        paths=[('baseline',ROOT/'results/matched_reconstruction_v6'/name/(name+'.xml'))]
        paths += [(v,ROOT/'models'/(name+'_matched_'+v+'_v6.xml')) for v in ['DNA_sequence_candidate','DNA_uniform_AA_surrogate_candidate']]
        for variant,p in paths:
            m=load(path=p);medium(m);bio=m.reactions.bio1
            for active,recycled in [('cpd11493_c0',['cpd12370_c0']),('cpd00166_c0',['cpd03422_c0','cpd01997_c0']),('cpd15665_c0',['cpd15666_c0'])]:
                assert all(abs(bio.get_coefficient(active)+bio.get_coefficient(k))<1e-9 for k in recycled)
                with m:
                    r=cobra.Reaction('V6_PAIRED_BIOMASS_MOIETY');r.bounds=(0,1000)
                    r.add_metabolites({m.metabolites.get_by_id(active):-1,**{m.metabolites.get_by_id(k):1 for k in recycled}});m.add_reactions([r]);m.objective=r
                    rows.append({'model':name,'variant':variant,'active_precursor':active,'recycled_pool':';'.join(recycled),'max_paired_demand':m.slim_optimize(error_value=np.nan),'status':m.solver.status,
                                 'meaning':'Produces the incorporated moiety with the same counterproduct as biomass; does not drain conserved R/Co carrier pool'})
                ledger.append({'model':name,'variant':variant,'active':active,'active_formula':m.metabolites.get_by_id(active).formula,'coefficient':bio.get_coefficient(active),'recycled':';'.join(recycled),'recycled_formula':';'.join(m.metabolites.get_by_id(k).formula for k in recycled),'recycled_coefficient':json.dumps({k:bio.get_coefficient(k) for k in recycled})})
    out=ROOT/'results/matched_biomass_v6';pd.DataFrame(rows).to_csv(out/'paired_carrier_precursor_demands.csv',index=False);pd.DataFrame(ledger).to_csv(out/'carrier_recycling_ledger.csv',index=False)
    record('results/matched_biomass_v6/precursor_test_interpretation.json',{'independent_demand_zero':['cpd11493_c0','cpd00166_c0','cpd15665_c0'],
        'reason':'The submitted biomass consumes active forms and returns apo-ACP, cobinamide plus dimethylbenzimidazole, and n-1 peptidoglycan in equal molar coefficients. Single-metabolite drains remove conserved generic R, cobalt or lower-ligand pools, unlike biomass. Paired demands test the incorporated chemical moiety while preserving the explicit template recycling convention. An initial incomplete cobalamin counterproduct test omitted dimethylbenzimidazole and returned zero; it is not accepted as the complete carrier-aware test.',
        'not_claimed':'This is a stoichiometric precursor-feasibility check, not measured cofactor composition, polymer chain length, carrier abundance or experimental biomass validation.'})
    print(pd.DataFrame(rows)[['model','variant','active_precursor','max_paired_demand']].to_string(index=False))
if __name__=='__main__':main()
