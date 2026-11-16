"""Native exchange FVA conditioned on the NEW certified native-QP growth rates."""
from native_micom_complete import *

def main():
    rates=pd.read_csv(OUT/'native_tradeoff_controls.csv');rows=[]
    for a in [.1,.25,.5,.75,.9]:
        c=community(a);c.solver='glpk';c.solver.configuration.tolerances.feasibility=1e-8
        for ctrl in CONTROLS:
            with c:
                control(c,ctrl)
                point=rates[(rates.M_abundance==a)&(rates.control==ctrl)&(rates.tradeoff==.95)].iloc[0]
                # Native QP rates before pFBA growth-preservation tolerance.
                p=pd.read_csv(OUT/f'qp_variables_a{a}_{ctrl}_f0.95_attempt0.csv',index_col=0).iloc[:,0].to_dict()
                for taxon in ['M','B']:
                    const=c.constraints['objective_'+taxon]
                    value=sum(float(coef)*(p[var.name] if hasattr(var,'name') else float(var)) for var,coef in const.expression.as_coefficients_dict().items())
                    const.ub=None;const.lb=max(0,value-1e-7);const.ub=value+1e-7
                ids=[f'EX_{code}_e0__{taxon}' for code in ['cpd00159','cpd00221','cpd00211','cpd00029','cpd00141'] for taxon in ['M','B'] if f'EX_{code}_e0__{taxon}' in c.reactions]
                c.objective=0
                f=flux_variability_analysis(c,reaction_list=ids,fraction_of_optimum=0,processes=1).reset_index(names='reaction')
                if f.isna().any().any():raise RuntimeError('Native conditioned FVA contains an unresolved optimum')
                f['M_abundance']=a;f['control']=ctrl;f['tradeoff']=.95;f['conditioning']='new certified native QP member growth, +/-1e-7; zero LP objective, no secondary community optimum constraint'
                rows.append(f)
                for taxon in ['M','B']:
                    const=c.constraints['objective_'+taxon];const.ub=None;const.lb=0
            pd.concat(rows).to_csv(OUT/'native_QP_conditioned_exchange_FVA.csv',index=False)
        print('native QP-conditioned FVA abundance',a,'complete',flush=True)

if __name__=='__main__':main()
