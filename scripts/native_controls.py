from native_analysis import *
def main():
    rows=[]
    for label in FILES:
        m=load(label);medium(m)
        for a in [.1,.25,.5,.75,.9]:
            for no_lactate in [False,True]:
                for allocation,factor in [('same_total_biomass_supply',1.),('all_pool_per_member_biomass',1/a)]:
                    with m:
                        for rid,cap in m.medium.items():m.reactions.get_by_id(rid).lower_bound=-cap*factor
                        if no_lactate:
                            for code in ['cpd00159','cpd00221']:
                                rid='EX_'+code+'_e0'
                                if rid in m.reactions:m.reactions.get_by_id(rid).lower_bound=0
                        rows.append({'model':label,'member_abundance':a,'no_external_lactate':no_lactate,'allocation':allocation,**phenotype(m)})
    pd.DataFrame(rows).to_csv(ROOT/'results/native_matched_monocultures.csv',index=False)
    probes={
        'ATP_hydrolysis':{'cpd00002_c0':-1,'cpd00001_c0':-1,'cpd00008_c0':1,'cpd00009_c0':1,'cpd00067_c0':1},
        'GTP_hydrolysis':{'cpd00038_c0':-1,'cpd00001_c0':-1,'cpd00031_c0':1,'cpd00009_c0':1,'cpd00067_c0':1},
        'UTP_hydrolysis':{'cpd00062_c0':-1,'cpd00001_c0':-1,'cpd00014_c0':1,'cpd00009_c0':1,'cpd00067_c0':1},
        'CTP_hydrolysis':{'cpd00052_c0':-1,'cpd00001_c0':-1,'cpd00096_c0':1,'cpd00009_c0':1,'cpd00067_c0':1},
        'NADH_oxidation':{'cpd00004_c0':-1,'cpd00003_c0':1,'cpd00067_c0':1},
        'NADPH_oxidation':{'cpd00005_c0':-1,'cpd00006_c0':1,'cpd00067_c0':1},
        'proton_influx':{'cpd00067_e0':-1,'cpd00067_c0':1},
        'sodium_influx':{'cpd00971_e0':-1,'cpd00971_c0':1}}
    rows=[]
    for p in sorted((ROOT/'models').glob('*.xml')):
        m=load(path=p)
        for r in m.boundary:r.bounds=(0,0)
        for name,stoich in probes.items():
            with m:
                absent=[x for x in stoich if x not in m.metabolites]
                if absent:rows.append({'model':p.stem,'probe':name,'status':'metabolites_missing','missing':';'.join(absent)});continue
                r=cobra.Reaction('V3_ENERGY_PROBE',lower_bound=0,upper_bound=1000);r.add_metabolites({m.metabolites.get_by_id(k):v for k,v in stoich.items()});m.add_reactions([r]);m.objective=r
                maximum=m.slim_optimize(error_value=np.nan);rows.append({'model':p.stem,'probe':name,'maximum':maximum,'status':m.solver.status,'stoichiometry':r.reaction})
    pd.DataFrame(rows).to_csv(ROOT/'results/native_closed_boundary_energy_probes.csv',index=False)
if __name__=='__main__':main()
