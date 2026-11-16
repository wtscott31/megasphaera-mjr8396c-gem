"""Promote sequence-supported annotation only; preserve reaction chemistry/GPRs."""
from native_analysis import *
from BOFdat.core import dna

def main():
    out=ROOT/'results/JGI_mapping';d=pd.read_csv(out/'deposited_model_477_gene_mapping.csv',dtype={'JGI_gene':str}).set_index('JGI_gene')
    m=load();before=mass(m.reactions.bio1.metabolites)
    for g in m.genes:
        row=d.loc[g.id]
        assert ';' not in row.public_protein_accessions
        g.annotation['ncbiprotein']=row.public_protein_accessions
        g.notes['public_locus_tag']=row.public_locus_tags
        g.notes['sequence_mapping_method']=row.mapping_status
        g.notes['original_protein_sha256']=row.AA_sha256
        g.notes['genomic_translation_validation']='Original protein translated exactly in original IMG contig and public assembly, allowing verified bacterial initiation codons'
    m.notes['revision_v4']='Sequence-supported gene cross-references only; reaction bounds, chemistry and GPR logic unchanged'
    dest=ROOT/'models/MJR8396C_sequence_annotated_primary.xml';cobra.io.write_sbml_model(m,str(dest))
    reread=load(path=dest);assert len(reread.genes)==477
    original=load()
    for r in original.reactions:
        rr=reread.reactions.get_by_id(r.id)
        assert rr.bounds==r.bounds and rr.gene_reaction_rule==r.gene_reaction_rule
        assert {x.id:v for x,v in rr.metabolites.items()}=={x.id:v for x,v in r.metabolites.items()}
    assert all(g.annotation.get('ncbiprotein') for g in reread.genes)
    medium(m);medium(original);assert abs(m.slim_optimize()-original.slim_optimize())<1e-10
    # Repeat the unmodified official DNA module on the actual original genome.
    bio=m.reactions.bio1;old={x.id:v for x,v in bio.metabolites.items()};ppi=m.metabolites.cpd00012_c0
    DNA_mass=sum(-old[DNA[x]]*(m.metabolites.get_by_id(DNA[x]).formula_weight-ppi.formula_weight) for x in DNA)/1000
    for source,target in MAP.items():m.metabolites.get_by_id(source).id=target
    m.repair();mapped=ROOT/'inputs/BOFdat_original_JGI_mapped_model.json';cobra.io.save_json_model(m,str(mapped))
    fasta=ROOT/'inputs/original_JGI_gapseq/Megasphaera_genome_assembly/2775506928.fna'
    raw=dna.generate_coefficients(str(fasta),str(mapped),DNA_WEIGHT_FRACTION=DNA_mass)
    raw={met.id:float(v) for met,v in raw.items()}
    rawmass=-sum(m.metabolites.get_by_id(k).formula_weight*v for k,v in raw.items())/1000
    letters=collections.Counter(''.join(str(r.seq).upper() for r in SeqIO.parse(fasta,'fasta')));n=sum(letters[x] for x in 'ACGT')
    f={'A':(letters['A']+letters['T'])/(2*n),'T':(letters['A']+letters['T'])/(2*n),'G':(letters['G']+letters['C'])/(2*n),'C':(letters['G']+letters['C'])/(2*n)}
    av=sum(f[x]*(m.metabolites.get_by_id(MAP[DNA[x]]).formula_weight-m.metabolites.ppi_c.formula_weight) for x in DNA)
    total=1000*DNA_mass/av;corrected={MAP[DNA[x]]:-total*f[x] for x in DNA};corrected['ppi_c']=total
    candidate=load(path=dest);b=candidate.reactions.bio1
    reverse={v:k for k,v in MAP.items()}
    changes={candidate.metabolites.get_by_id(DNA[x]):corrected[MAP[DNA[x]]]-old[DNA[x]] for x in DNA}
    changes[candidate.metabolites.cpd00012_c0]=total-sum(-old[DNA[x]] for x in DNA)
    b.add_metabolites(changes)
    assert abs(mass(b.metabolites)-before)<1e-10
    cobra.io.write_sbml_model(candidate,str(ROOT/'models/MJR8396C_original_JGI_DNA_mass_corrected_candidate.xml'))
    record('results/JGI_mapping/annotation_and_original_DNA_audit.json',{'annotation_only_primary':str(dest.relative_to(ROOT)),'all_477_genes_cross_referenced':True,'reaction_chemistry_bounds_GPR_unchanged_roundtrip_verified':True,'original_JGI_DNA_weight_fraction_assumed':DNA_mass,'native_official_BOFdat_raw_coefficients':raw,'native_raw_DNA_mass':rawmass,'native_fraction_recovered':rawmass/DNA_mass,'native_raw_rejected':True,'external_mass_corrected_coefficients':corrected,'genome_derived_double_strand_frequencies':f,'candidate_net_mass':mass(b.metabolites),'baseline_net_mass':before,'native_source_sha256':sha(Path(dna.__file__))})
    # Connect all existing locus-level functional evidence to actual model IDs.
    allmap=pd.read_csv(out/'all_original_protein_mapping.csv',dtype={'JGI_gene':str})
    links=allmap.assign(public_protein_accessions=allmap.public_protein_accessions.fillna('').str.split(';')).explode('public_protein_accessions')
    for name in ['alignment_evidence.csv','candidate_neighborhoods_10kb.csv','complex_component_candidate_inventory.csv']:
        src=ROOT/'results/locus_evidence'/name;df=pd.read_csv(src)
        cols=[c for c in df.columns if 'protein' in c.lower() or 'subject' in c.lower()]
        print(name,cols,flush=True)
    record('results/JGI_mapping/promoted_edits.json',{'promoted':'477 gene cross-references and traceable locus notes','not_promoted':'New enzyme activity, reaction directions, stoichiometry or GPR logic','DNA':'Separate original-genome sequence/mass-corrected sensitivity candidate; original presumed dry-weight fraction retained'})
    print('Annotation-only primary and original-JGI DNA candidate validated',flush=True)

if __name__=='__main__':main()
