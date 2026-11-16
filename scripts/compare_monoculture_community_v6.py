"""Direct effective-supply mono/community comparison with matched-growth FVA."""
from matched_micom_v6 import pairs
from native_analysis import *

def main():
    rows=[]
    for pair,(mf,bf) in pairs().items():
        d=ROOT/'results/matched_micom_v6'/pair
        co=pd.read_csv(d/'native_tradeoff_controls.csv');fl=pd.read_csv(d/'native_taxon_community_fluxes.csv')
        m=load(path=mf)
        for a in [.1,.25,.5,.75,.9]:
            for ctrl in ['baseline','no_external_lactate','no_lactate_block_both']:
                selected=co[co.M_abundance.eq(a)&co.tradeoff.eq(.95)&co.control.eq(ctrl)].iloc[0]
                cf=fl[fl.M_abundance.eq(a)&fl.tradeoff.eq(.95)&fl.control.eq(ctrl)&fl.compound.eq('cpd00211')].iloc[0]
                with m:
                    medium(m);m.medium={k:v/a for k,v in m.medium.items()}
                    if ctrl!='baseline':
                        for code in ['cpd00159','cpd00221']:
                            rid='EX_'+code+'_e0'
                            if rid in m.reactions:m.reactions.get_by_id(rid).lower_bound=0
                    mono=phenotype(m);r=next(iter(cobra.util.solver.linear_reaction_coefficients(m)))
                    lo=hi=np.nan;status='unmatched_growth';g=float(selected.growth_M)
                    if mono['growth']+1e-6>=g and 'EX_cpd00211_e0' in m.reactions:
                        r.bounds=(max(0,g-1e-7),g+1e-7)
                        f=flux_variability_analysis(m,reaction_list=['EX_cpd00211_e0'],fraction_of_optimum=0,processes=1)
                        lo,hi=map(float,f.loc['EX_cpd00211_e0',['minimum','maximum']]);status='optimal'
                    rows.append({'pair':pair,'M_abundance':a,'tradeoff':.95,'control':ctrl,'monoculture_growth_max':mono['growth'],
                                 'community_member_growth':g,'growth_difference_co_minus_mono':g-mono['growth'],'monoculture_pFBA_butyrate_at_max_growth':mono.get('cpd00211',np.nan),
                                 'community_member_butyrate':cf.M_specific,'mono_butyrate_min_at_matched_member_growth':lo,'mono_butyrate_max_at_matched_member_growth':hi,
                                 'matched_growth_FVA_status':status,'co_butyrate_inside_mono_interval':bool(lo-1e-5<=cf.M_specific<=hi+1e-5) if np.isfinite(hi) else None,
                                 'supply_basis':'Same pooled nutrient capacity divided by M abundance; all organic imports from the stated medium. Partner-derived metabolites unavailable to monoculture.'})
        print(pair,'effective-supply comparison complete',flush=True)
    f=pd.DataFrame(rows);f.to_csv(ROOT/'tables/S80b_effective_supply_mono_community_v6.csv',index=False)
    f[f.M_abundance.eq(.5)].to_csv(ROOT/'results/reviewer_v6/effective_supply_selected_comparison.csv',index=False)
if __name__=='__main__':main()
