"""Native COBRApy reruns and DNA-only candidate with explicit mass correction."""
from common import *
import os,re,collections,warnings
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'logs/matplotlib'))
import cobra,numpy as np,pandas as pd
from Bio import SeqIO
from cobra.flux_analysis import pfba,flux_variability_analysis
cobra.Configuration().solver='glpk'
cobra.Configuration().processes=1
FILES={'M':BASE/'models/Megasphaera_sp_MJR8396C_filled_anaerobic.xml','B':BASE/'models/Bifidobacterium_adolescentis_filled_anaerobic.xml'}
DNA={'A':'cpd00115_c0','T':'cpd00357_c0','G':'cpd00241_c0','C':'cpd00356_c0'}
MAP={DNA[x]:{'A':'datp_c','T':'dttp_c','G':'dgtp_c','C':'dctp_c'}[x] for x in DNA}|{'cpd00012_c0':'ppi_c'}
MW={'H':1.00794,'C':12.0107,'N':14.0067,'O':15.9994,'P':30.973762,'S':32.065,'Fe':55.845,'K':39.0983,'Ca':40.078,'Mn':54.938045,'Zn':65.38,'Co':58.933195,'Cu':63.546,'Mg':24.305,'Cl':35.453,'Na':22.989769}
TARGETS=['cpd00027','cpd00159','cpd00221','cpd00029','cpd00211','cpd00141','cpd00047','cpd00036','cpd00013','cpd00239','cpd00324','cpd00011','cpd11640']
def load(label='M',path=None):
    m=cobra.io.read_sbml_model(str(path or FILES[label]));m.solver='glpk';m.solver.configuration.tolerances.feasibility=1e-8;return m
def medium(m,name='gut_like'):
    df=pd.read_csv(BASE/f'media/MJR8396C_{name}_medium.tsv',sep='\t');available={r.id:r for r in m.exchanges};mapped={};missing=[]
    for row in df.itertuples():
        rid='EX_'+row.compounds+'_e0'
        if rid in available:mapped[rid]=float(row.maxFlux)
        else:missing.append(row.compounds)
    m.medium={k:v for k,v in mapped.items() if v>0};return missing
def mass(stoich):
    atoms=collections.defaultdict(float)
    for met,v in stoich.items():
        if met.id=='cpd11416_c0':continue
        for element,n in met.elements.items():atoms[element]-=v*n
    unknown={k:v for k,v in atoms.items() if k not in MW and abs(v)>1e-8}
    if unknown:raise ValueError(unknown)
    return sum(MW.get(k,0)*v for k,v in atoms.items())/1000
def phenotype(m):
    s=m.optimize();row={'status':s.status,'growth':s.objective_value}
    if s.status=='optimal':
        p=pfba(m,fraction_of_optimum=1-1e-9)
        row.update({cpd:float(p.fluxes.get('EX_'+cpd+'_e0',np.nan)) for cpd in TARGETS})
        S=cobra.util.array.create_stoichiometric_matrix(m,array_type='lil').tocsr()
        row['mass_residual']=float(np.max(np.abs(S@p.fluxes.values)))
    return row
def baseline():
    rows=[];fvas=[]
    for label in FILES:
        m=load(label)
        for med in ['gut_like','carbon_rich','AA_rich']:
            missing=medium(m,med);rows.append({'model':label,'medium':med,'unmapped_medium':';'.join(missing),**phenotype(m)})
            ids=['EX_'+x+'_e0' for x in TARGETS if 'EX_'+x+'_e0' in m.reactions]
            for fraction in [.9,.99,1.]:
                f=flux_variability_analysis(m,reaction_list=ids,fraction_of_optimum=fraction,processes=1)
                f=f.reset_index(names='reaction');f['model']=label;f['medium']=med;f['growth_fraction']=fraction;fvas.append(f)
        pd.DataFrame(rows).to_csv(ROOT/'results/native_COBRA_baselines.csv',index=False)
        pd.concat(fvas).to_csv(ROOT/'results/native_COBRA_exchange_FVA.csv',index=False)
    old=pd.read_csv(BASE/'extended/results/strain_design/minimal_gene_sets_and_global_GPR_knockouts.csv').fillna('');rows=[];m=load();medium(m)
    for row in old.itertuples():
        with m:
            genes=str(row.genes).split(';') if row.genes else []
            cobra.manipulation.knock_out_model_genes(m,genes)
            g=m.slim_optimize(error_value=float('nan'));rows.append({'design':row.design,'genes':';'.join(genes),'native_growth':g,'v2_growth':row.growth,'difference':g-float(row.growth),'solver_status':m.solver.status})
    pd.DataFrame(rows).to_csv(ROOT/'results/native_COBRA_GPR_knockouts.csv',index=False)
    old=pd.read_csv(BASE/'extended/results/strain_design/historical_reaction_deletion_retests.csv').fillna('');rows=[]
    for row in old.itertuples():
        with m:
            ids=str(row.target_reactions).split(';') if row.target_reactions else []
            for rid in ids:m.reactions.get_by_id(rid).knock_out()
            g=m.slim_optimize(error_value=float('nan'));rows.append({'design':row.design,'targets':';'.join(ids),'native_growth':g,'v2_growth':row.growth,'difference':g-float(row.growth),'solver_status':m.solver.status})
    pd.DataFrame(rows).to_csv(ROOT/'results/native_COBRA_reaction_knockouts.csv',index=False)
def biomass():
    from BOFdat.core import dna
    m=load();b=m.reactions.bio1;old={x.id:v for x,v in b.metabolites.items()};ppi=m.metabolites.cpd00012_c0
    DNA_mass=sum(-old[DNA[x]]*(m.metabolites.get_by_id(DNA[x]).formula_weight-ppi.formula_weight) for x in DNA)/1000
    fasta=ROOT/'inputs/genomes/GCA_001546855.1/genomic.fna'
    for source,target in MAP.items():m.metabolites.get_by_id(source).id=target
    m.repair();dest=ROOT/'inputs/BOFdat_mapped_model.json';cobra.io.save_json_model(m,str(dest))
    native=dna.generate_coefficients(str(fasta),str(dest),DNA_WEIGHT_FRACTION=DNA_mass)
    raw={met.id:float(v) for met,v in native.items()};rawmass=-sum(m.metabolites.get_by_id(k).formula_weight*v for k,v in raw.items())/1000
    letters=collections.Counter(''.join(str(r.seq).upper() for r in SeqIO.parse(fasta,'fasta')));n=sum(letters[x] for x in 'ACGT')
    freq={'A':(letters['A']+letters['T'])/(2*n),'T':(letters['A']+letters['T'])/(2*n),'G':(letters['G']+letters['C'])/(2*n),'C':(letters['G']+letters['C'])/(2*n)}
    average=sum(freq[x]*(m.metabolites.get_by_id(MAP[DNA[x]]).formula_weight-m.metabolites.ppi_c.formula_weight) for x in DNA)
    total=1000*DNA_mass/average;corrected={MAP[DNA[x]]:-total*freq[x] for x in DNA};corrected['ppi_c']=total
    corrected_mass=-sum(m.metabolites.get_by_id(k).formula_weight*v for k,v in corrected.items())/1000
    assert abs(corrected_mass-DNA_mass)<1e-12
    assert abs(sum(corrected[k] for k in corrected if k!='ppi_c')+corrected['ppi_c'])<1e-12
    record('results/BOFdat_native_and_mass_correction.json',{'official_repository':'https://github.com/jclachance/BOFdat','commit':'cf2de074c2d789ad47b890083cbbffadab40f177','native_source_unmodified':True,'mapping':MAP,'DNA_fraction_g_gDW':DNA_mass,'fraction_provenance':'Preserved original deposited biomass DNA module mass; assumption, not a measured organism DNA fraction.','sequence_frequencies_double_strand':freq,'native_coefficients':raw,'native_net_mass_g_gDW':rawmass,'native_fraction_of_requested_mass':rawmass/DNA_mass,'correction_coefficients':corrected,'corrected_mass_g_gDW':corrected_mass,'correction_method':'Native _get_ratio returns ratios summing to 0.5 after double-strand counting. Native conversion additionally treats base count fractions as per-base mass fractions. Corrected coefficients use double-strand nucleotide mole fractions and one mean polymer residue mass (model dNTP minus exact model PPi), mmol/gDW = 1000 * assumed DNA mass fraction / mean residue g/mol. This transparent external correction is NOT represented as unmodified native BOFdat output.','limits':'DNA-only, model formula/protonation convention; negligible polymer end groups ignored; no measured RNA/proteomics/lipids/GAM/NGAM inferred.'})
    # Restore original identifiers before any exported model.
    for source,target in MAP.items():m.metabolites.get_by_id(target).id=source
    m.repair();cobra.io.write_sbml_model(m,str(ROOT/'models/MJR8396C_baseline.xml'))
    reverse={v:k for k,v in MAP.items()};rows=[];variants=[('baseline',None),('DNA_native_raw_diagnostic',raw),('DNA_sequence_mass_corrected_candidate',corrected)]
    summaries=[];prec=[];sens=[];energy=[]
    for name,coeff in variants:
        v=m.copy();vb=v.reactions.bio1
        if coeff:
            delta={v.metabolites.get_by_id(source):coeff[MAP[source]]-old[source] for source in DNA.values()}
            old_DNA_ppi=-sum(old[x] for x in DNA.values())
            delta[v.metabolites.cpd00012_c0]=coeff['ppi_c']-old_DNA_ppi;vb.add_metabolites(delta)
        if name!='baseline':cobra.io.write_sbml_model(v,str(ROOT/'models'/f'MJR8396C_{name}.xml'))
        for met,c in vb.metabolites.items():rows.append({'variant':name,'metabolite':met.id,'coefficient':c,'baseline_coefficient':old[met.id],'delta':c-old[met.id]})
        total_mass=mass(vb.metabolites)
        for med in ['gut_like','carbon_rich','AA_rich']:
            medium(v,med);summaries.append({'variant':name,'medium':med,'net_biomass_mass_g_per_unit':total_mass,**phenotype(v)})
        medium(v)
        for substrate in ['cpd00159','cpd00221']:
            for amino in [True,False]:
                with v:
                    for c in ['cpd00027','cpd00029','cpd00159','cpd00221']:
                        if 'EX_'+c+'_e0' in v.reactions:v.reactions.get_by_id('EX_'+c+'_e0').lower_bound=0
                    if not amino:
                        for c in ['cpd00041','cpd00033','cpd00039','cpd00060','cpd00084','cpd00119']:v.reactions.get_by_id('EX_'+c+'_e0').lower_bound=0
                    if 'EX_'+substrate+'_e0' in v.reactions:v.reactions.get_by_id('EX_'+substrate+'_e0').lower_bound=-10
                    summaries.append({'variant':name,'medium':'lactate_conditioned','substrate':substrate,'AA_background':amino,'net_biomass_mass_g_per_unit':total_mass,**phenotype(v)})
        for gam in [20,40,60]:
            for ngam in [0,1,2,5]:
                with v:
                    current=-vb.get_coefficient('cpd00002_c0') # Includes RNA ATP; GAM is explicitly 40 from template ADP.
                    change=gam-40
                    vb.add_metabolites({v.metabolites.get_by_id(k):change*s for k,s in {'cpd00002_c0':-1,'cpd00001_c0':-1,'cpd00008_c0':1,'cpd00009_c0':1,'cpd00067_c0':1}.items()})
                    v.reactions.rxn00062_c0.lower_bound=ngam
                    sens.append({'variant':name,'GAM_assumed':gam,'NGAM_assumed':ngam,**phenotype(v)})
        for met,c in list(vb.metabolites.items()):
            if c>=0:continue
            with v:
                demand=v.add_boundary(met,type='demand',reaction_id='V3_DEMAND',lb=0,ub=1000);v.objective=demand
                val=v.slim_optimize(error_value=np.nan);prec.append({'variant':name,'metabolite':met.id,'max_demand':val,'solver_status':v.solver.status})
        with v:
            for r in v.boundary:r.bounds=(0,0)
            v.objective=v.reactions.rxn00062_c0
            energy.append({'variant':name,'probe':'closed-boundary ATP hydrolysis','maximum':v.slim_optimize(error_value=np.nan),'status':v.solver.status})
        pd.DataFrame(rows).to_csv(ROOT/'results/biomass_coefficients.csv',index=False)
        pd.DataFrame(summaries).to_csv(ROOT/'results/biomass_candidate_phenotypes.csv',index=False)
        pd.DataFrame(prec).to_csv(ROOT/'results/biomass_precursor_tests.csv',index=False)
        pd.DataFrame(sens).to_csv(ROOT/'results/biomass_maintenance_sensitivity.csv',index=False)
        pd.DataFrame(energy).to_csv(ROOT/'results/biomass_ATP_cycle_tests.csv',index=False)
    # Genome-derived uniform-protein surrogate; no biomass protein replacement.
    counts=collections.Counter(''.join(str(r.seq) for r in SeqIO.parse(ROOT/'inputs/genomes/GCA_001546855.1/protein.faa','fasta')))
    pd.DataFrame([{'AA':x,'count':counts[x],'frequency':counts[x]/sum(counts.values()),'provenance':'Uniform abundance across annotated CDS; genome-derived surrogate, not proteomics'} for x in sorted(counts)]).to_csv(ROOT/'results/protein_uniform_sequence_surrogate.csv',index=False)
if __name__=='__main__':
    baseline();biomass()
