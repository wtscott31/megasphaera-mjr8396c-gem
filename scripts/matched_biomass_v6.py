"""Sequence-informed biomass sensitivity; no target measurements invented."""
from native_analysis import *
from empirical_biomass_evidence import AA
import gzip,itertools

def main():
    out=ROOT/'results/matched_biomass_v6';out.mkdir(exist_ok=True);rows=[];ledger=[];prec=[];energy=[];assumptions=[]
    for n in ['MJR8396C','NM10','BL7']:
        d=ROOT/'results/matched_reconstruction_v6'/n;m=load(path=d/(n+'.xml'));bio=m.reactions.bio1;old={x.id:v for x,v in bio.metabolites.items()}
        oldmass=mass(bio.metabolites);ppi=m.metabolites.cpd00012_c0;water=m.metabolites.cpd00001_c0
        dm=sum(-old[DNA[b]]*(m.metabolites.get_by_id(DNA[b]).formula_weight-ppi.formula_weight) for b in DNA)/1000
        counts=collections.Counter(''.join(str(x.seq).upper() for x in SeqIO.parse(d/(n+'.fna'),'fasta')));den=sum(counts[b] for b in 'ATGC');gc=(counts['G']+counts['C'])/den
        freq={'A':(1-gc)/2,'T':(1-gc)/2,'G':gc/2,'C':gc/2}
        mean=sum(freq[b]*(m.metabolites.get_by_id(DNA[b]).formula_weight-ppi.formula_weight) for b in DNA)
        new={DNA[b]:-1000*dm*freq[b]/mean for b in DNA}
        changes={m.metabolites.get_by_id(k):v-old[k] for k,v in new.items()}
        changes[ppi]=-sum(new.values())+sum(old[DNA[b]] for b in DNA)
        bio.add_metabolites(changes);assert abs(mass(bio.metabolites)-oldmass)<1e-9
        dna=m.copy();protein=sum(-old[c+'_c0']*(m.metabolites.get_by_id(c+'_c0').formula_weight-water.formula_weight) for c in AA.values())/1000
        with gzip.open(d/(n+'.faa.gz'),'rt') as f:ac=collections.Counter(''.join(str(x.seq).rstrip('*') for x in SeqIO.parse(f,'fasta')))
        af={a:ac[a] for a in AA};mean_aa=sum(v*(m.metabolites.get_by_id(AA[a]+'_c0').formula_weight-water.formula_weight) for a,v in af.items())
        sto={AA[a]+'_c0':-1000*protein*v/mean_aa for a,v in af.items()}
        before_n=sum(-bio.get_coefficient(c+'_c0') for c in AA.values());after_n=-sum(sto.values())
        changes={m.metabolites.get_by_id(k):v-bio.get_coefficient(k) for k,v in sto.items()};changes[water]=after_n-before_n
        bio.add_metabolites(changes);assert abs(mass(bio.metabolites)-oldmass)<1e-9
        assumptions.append({'model':n,'GC_fraction':gc,'assumed_DNA_mass':dm,'assumed_protein_mass':protein,'template_total_mass':oldmass,'DNA_evidence':'Genome-derived double-strand mole fractions; inherited DNA mass fraction','protein_evidence':'Uniform abundance across translated ORFs; NOT measured proteomics','remaining':'RNA/lipid/cell-wall fractions and GAM/NGAM remain template assumptions'})
        for label,v in [('DNA_sequence_candidate',dna),('DNA_uniform_AA_surrogate_candidate',m)]:
            v.id=n+'_matched_'+label+'_v6';dest=ROOT/'models'/(v.id+'.xml');cobra.io.write_sbml_model(v,str(dest))
            b=v.reactions.bio1
            for met,c in b.metabolites.items():ledger.append({'model':n,'variant':label,'metabolite':met.id,'coefficient':c,'baseline':old[met.id]})
            for med,gam,ngam in itertools.product(['gut_like','carbon_rich','AA_rich'],[20,40,60],[0,1,2,5]):
                with v:
                    medium(v,med);delta=gam-b.get_coefficient('cpd00008_c0')
                    b.add_metabolites({v.metabolites.get_by_id(c+'_c0'):delta*k for c,k in {'cpd00002':-1,'cpd00001':-1,'cpd00008':1,'cpd00009':1,'cpd00067':1}.items()})
                    v.reactions.rxn00062_c0.lower_bound=ngam
                    rows.append({'model':n,'variant':label,'medium':med,'GAM_assumed':gam,'NGAM_assumed':ngam,'biomass_mass':mass(b.metabolites),**phenotype(v)})
            medium(v)
            for met,c in list(b.metabolites.items()):
                if c>=0:continue
                with v:
                    dmrx=v.add_boundary(met,type='demand',reaction_id='V6_BIOMASS_PRECURSOR');v.objective=dmrx
                    prec.append({'model':n,'variant':label,'metabolite':met.id,'max_demand':v.slim_optimize(error_value=np.nan),'status':v.solver.status})
            with v:
                for r in v.boundary:r.bounds=(0,0)
                v.reactions.rxn00062_c0.lower_bound=0;v.objective=v.reactions.rxn00062_c0
                energy.append({'model':n,'variant':label,'ATP_max':v.slim_optimize(error_value=np.nan),'status':v.solver.status})
        print(n,'sequence-biomass sensitivities complete',flush=True)
    for file,data in [('growth_maintenance_sensitivity',rows),('coefficients',ledger),('precursor_producibility',prec),('ATP_cycle_checks',energy)]:pd.DataFrame(data).to_csv(out/(file+'.csv'),index=False)
    record('results/matched_biomass_v6/assumptions.json',assumptions)
if __name__=='__main__':main()
