"""Trace primary related-strain measurements and test a labelled proxy BOF.

No measurement is attributed to MJR8396C. LC1 hydrolysate amino-acid ratios
inform a sensitivity model's shape; the inherited total protein mass and
unmeasured Cys/Trp remain assumptions. Polymerization water is mass balanced.
"""
from native_analysis import *
import zipfile,lxml.etree as E

AA={'D':'cpd00041','N':'cpd00132','E':'cpd00023','Q':'cpd00053','M':'cpd00060','S':'cpd00054','A':'cpd00035','G':'cpd00033','K':'cpd00039','R':'cpd00051','W':'cpd00065','F':'cpd00066','Y':'cpd00069','C':'cpd00084','T':'cpd00161','L':'cpd00107','H':'cpd00119','V':'cpd00156','P':'cpd00129','I':'cpd00322'}
SOURCE='https://doi.org/10.1128/aem.51.5.1141-1143.1986'

def main():
    out=ROOT/'results/empirical_biomass';out.mkdir(exist_ok=True)
    # Direct transcription, visually checked against Wallace 1986 Table 1,
    # p.1142; column "M. elsdenii cells", micromol/g dry weight, means of 3.
    measured={'Asx':420,'T':217,'S':201,'Glx':470,'G':380,'A':488,'M':109,'V':245,'I':171,'L':277,'Y':114,'F':141,'H':86,'K':254,'R':158,'P':188}
    pd.DataFrame([{'pool':k,'micromol_per_gDW':v,'organism':'Megasphaera elsdenii LC1','strain_identifier':'NCIB 8927','determinations':3,'source':SOURCE,'source_table':'Table 1, p.1142','evidence':'measured whole-cell acid hydrolysate; related strain, not MJR8396C'} for k,v in measured.items()]).to_csv(out/'Wallace1986_measured_amino_acid_pools.csv',index=False)
    lipids=[(27,71.4,28.1,46,None,None),(30,68.4,31.5,57,1.8,3),(37,65.2,34.8,49.5,3,4),(42.5,64.7,35.3,54,None,2)]
    pd.DataFrame([{'temperature_C':t,'PE_percent_lipid_P':pe,'PS_percent_lipid_P':ps,'plasmalogen_alkenyl_bonds_per_lipid_P_percent':pl,'PE_SE':se,'plasmalogen_SE':pse,'organism':'Megasphaera elsdenii ATCC 17752','source':'https://doi.org/10.1128/jb.149.2.567-575.1982','table':'Table 4, p.571','evidence':'related-strain lipid composition, not whole-cell lipid dry-mass fraction'} for t,pe,ps,pl,se,pse in lipids]).to_csv(out/'Johnston1982_measured_lipid_headgroups.csv',index=False)
    # Maintenance source Table 1, p.539, with an ammonia correction that itself
    # assumes carbohydrate equivalence of amino-acid catabolism.
    D=[.410,.319,.223,.196,.154,.094];total=[.142,.244,.330,.307,.326,.357];dry=[.053,.099,.104,.098,.094,.087]
    maint=pd.DataFrame({'dilution_h_inverse':D,'carbohydrate_equivalent_mg_ml':total,'dry_weight_mg_ml':dry})
    maint['yield_gDW_g_carbohydrate']=maint.dry_weight_mg_ml/maint.carbohydrate_equivalent_mg_ml
    fit=np.polyfit(1/maint.dilution_h_inverse,1/maint.yield_gDW_g_carbohydrate,1)
    maint['source']='https://doi.org/10.1128/aem.37.3.537-543.1979';maint['strain']='M. elsdenii B159';maint.to_csv(out/'Russell1979_measured_maintenance_table.csv',index=False)
    record('results/empirical_biomass/maintenance_fit.json',{'method':'ordinary least squares of inverse yield against inverse dilution, rounded table entries','slope_g_carbohydrate_per_gDW_h':fit[0],'intercept_inverse_yield':fit[1],'max_yield_gDW_per_g_carbohydrate':1/fit[1],'published_slope':.187,'published_max_yield_approx':.42,'ATP_NGAM_not_inferred':'Requires ATP yield, amino-acid energy assumptions, strain and culture transfer; not a measured MJR8396C ATP drain'})
    p=ROOT/'inputs/comparators/iME375/Supplementary data 1.docx'
    with zipfile.ZipFile(p) as z:
        document=E.fromstring(z.read('word/document.xml'))
        paragraphs=[''.join(el.xpath('.//*[local-name()="t"]/text()')) for el in document.xpath('//*[local-name()="p"]')]
    (out/'iME375_S1_extracted_text.txt').write_text('\n'.join(paragraphs))
    m=load(path=ROOT/'models/MJR8396C_original_JGI_DNA_mass_corrected_candidate.xml');b=m.reactions.bio1
    old={x.id:v for x,v in b.metabolites.items()};water=m.metabolites.cpd00001_c0
    protein=sum(-old[mid+'_c0']*(m.metabolites.get_by_id(mid+'_c0').formula_weight-water.formula_weight) for mid in AA.values())/1000
    fixed_protein=sum(-old[AA[aa]+'_c0']*(m.metabolites.get_by_id(AA[aa]+'_c0').formula_weight-water.formula_weight) for aa in ['C','W'])/1000
    relative={k:v/1000 for k,v in measured.items() if k not in ['Asx','Glx']}
    relative.update(D=.210,N=.210,E=.235,Q=.235)
    shape_mass=sum(relative[a]*(m.metabolites.get_by_id(AA[a]+'_c0').formula_weight-water.formula_weight) for a in relative)/1000
    scale=(protein-fixed_protein)/shape_mass
    changes={m.metabolites.get_by_id(AA[a]+'_c0'):-value*scale-old[AA[a]+'_c0'] for a,value in relative.items()}
    new_n=sum(relative.values())*scale+sum(-old[AA[a]+'_c0'] for a in ['C','W'])
    old_n=sum(-old[mid+'_c0'] for mid in AA.values());changes[water]=new_n-old_n
    before=mass(b.metabolites);b.add_metabolites(changes);after=mass(b.metabolites)
    assert abs(after-before)<1e-10
    m.id='MJR8396C_related_strain_AA_shape_DNA_candidate'
    m.notes['evidence_scope']='Sensitivity candidate: LC1 measured amino-acid shape; assumed protein fraction, Asx/Glx 50:50 split, Cys/Trp inherited; not empirical MJR8396C biomass'
    cobra.io.write_sbml_model(m,str(ROOT/'models/MJR8396C_related_strain_AA_shape_DNA_candidate.xml'))
    record('results/empirical_biomass/candidate_assumptions.json',{'primary_model_unchanged':True,'candidate_only':True,'related_strain':'M. elsdenii LC1 / NCIB 8927','target_strain_measurements_available':False,'measured_input':'16 whole-cell amino-acid pools; original proportions, not iME375 amino-acid mass conversions','protein_mass_assumed_g_per_gDW':protein,'Cys_Trp':'inherited original coefficients; unmeasured','Asx_Glx':'50:50 molar split into acid/amide, inherited assumption','shape_scaling':scale,'before_total_mass':before,'after_total_mass':after,'RNA_lipid_GAM_NGAM':'unchanged assumptions','lipid_not_replaced':'Measured related-strain PS/plasmalogens are not represented equivalently in the original nine PG/CL/PE species; no speculative pathway addition','release_gate':'Does not resolve missing strain-specific whole-cell fractions or maintenance measurements'})
    pd.DataFrame([{'metabolite':x.id,'name':x.name,'baseline_DNA_candidate':old[x.id],'empirical_shape_candidate':v,'delta':v-old[x.id]} for x,v in b.metabolites.items()]).to_csv(out/'candidate_coefficients.csv',index=False)
    rows=[];fv=[];demands=[];energy=[]
    for name,path in [('original',FILES['M']),('original_JGI_DNA_candidate',ROOT/'models/MJR8396C_original_JGI_DNA_mass_corrected_candidate.xml'),('related_strain_AA_shape_DNA_candidate',ROOT/'models/MJR8396C_related_strain_AA_shape_DNA_candidate.xml')]:
        v=load(path=path)
        for med in ['gut_like','carbon_rich','AA_rich']:
            medium(v,med);rows.append({'model':name,'medium':med,**phenotype(v)})
            f=flux_variability_analysis(v,reaction_list=['EX_'+c+'_e0' for c in ['cpd00211','cpd00141','cpd00029','cpd00159','cpd00221']],fraction_of_optimum=.95,processes=1).reset_index(names='reaction');f['model']=name;f['medium']=med;fv.append(f)
        medium(v)
        for substrate in ['cpd00159','cpd00221']:
            for background in [True,False]:
                with v:
                    for code in ['cpd00027','cpd00029','cpd00159','cpd00221']:v.reactions.get_by_id('EX_'+code+'_e0').lower_bound=0
                    if not background:
                        for code in ['cpd00041','cpd00033','cpd00039','cpd00060','cpd00084','cpd00119']:v.reactions.get_by_id('EX_'+code+'_e0').lower_bound=0
                    v.reactions.get_by_id('EX_'+substrate+'_e0').lower_bound=-10
                    rows.append({'model':name,'medium':'lactate_conditioned','substrate':substrate,'AA_background':background,**phenotype(v)})
        for met,coef in list(v.reactions.bio1.metabolites.items()):
            if coef>=0:continue
            with v:
                dm=v.add_boundary(met,type='demand');v.objective=dm;value=v.slim_optimize(error_value=np.nan)
                demands.append({'model':name,'metabolite':met.id,'demand_max':value,'status':v.solver.status})
        for med in ['gut_like','carbon_rich','AA_rich']:
            medium(v,med)
            for gam in [20,40,60]:
                for ngam in [0,1,2,5]:
                    with v:
                        delta=gam-40
                        v.reactions.bio1.add_metabolites({v.metabolites.cpd00002_c0:-delta,v.metabolites.cpd00001_c0:-delta,v.metabolites.cpd00008_c0:delta,v.metabolites.cpd00009_c0:delta,v.metabolites.cpd00067_c0:delta})
                        atp=v.reactions.get_by_id('rxn00062_c0') if 'rxn00062_c0' in v.reactions else None
                        if atp is None:raise ValueError('Missing maintenance reaction')
                        atp.lower_bound=ngam
                        rows.append({'model':name,'medium':med,'GAM':gam,'NGAM':ngam,**phenotype(v)})
        for bounds in ['all_closed','inorganics_only']:
            with v:
                v.medium={}
                if bounds=='inorganics_only':v.medium={r.id:100 for r in v.exchanges if next(iter(r.metabolites)).id.split('_')[0] in ['cpd00001','cpd00009','cpd00011','cpd00013','cpd00067','cpd00048']}
                atp=v.reactions.get_by_id('rxn00062_c0');atp.lower_bound=0;v.objective=atp
                energy.append({'model':name,'condition':bounds,'ATP_max':v.slim_optimize(),'status':v.solver.status})
    pd.DataFrame(rows).to_csv(out/'phenotype_and_maintenance_sensitivity.csv',index=False)
    pd.concat(fv).to_csv(out/'selected_exchange_FVA.csv',index=False)
    pd.DataFrame(demands).to_csv(out/'precursor_demands.csv',index=False)
    pd.DataFrame(energy).to_csv(out/'ATP_cycle_checks.csv',index=False)
    print('Empirical evidence and sensitivity models executed',len(rows),'growth/maintenance rows',flush=True)

if __name__=='__main__':main()
