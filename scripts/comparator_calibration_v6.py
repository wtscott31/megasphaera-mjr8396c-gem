"""Bounded comparator calibration: common rich medium + physical mass audits.

No native comparator is gapfilled. New candidate repairs are separate models.
Unknown polymer mass is retained as unknown, never assigned a convenient mass.
"""
from comparators import *
from empirical_biomass_evidence import AA

OUT=ROOT/'results/comparator_calibration_v6'
AA_VMH=['asp_L','asn_L','glu_L','gln_L','met_L','ser_L','ala_L','gly','lys_L','arg_L','trp_L','phe_L','tyr_L','cys_L','thr_L','leu_L','his_L','val_L','pro_L','ile_L']
AA_SEED=dict(zip(AA_VMH,AA.values()))
IME_AA=dict(zip(['gly','ala_L','val_L','leu_L','ile_L','glu_L','gln_L','met_L','cys_L','asp_L','asn_L','pro_L','trp_L','tyr_L','his_L','phe_L','ser_L','thr_L','lys_L','arg_L'],['R'+str(i).zfill(3) for i in range(22,42)]))

def files():
    src=json.loads((ROOT/'results/comparators/structural_native_settings_audit.json').read_text())
    return {x['model']:ROOT/x['file'] for x in src}

def ime_candidate():
    m=load(path=ROOT/'inputs/comparators/iME375/Supplementary data 4.xml');ledger=[]
    # Repair only modules for which chemical monomers and source ratios are
    # explicit. Lipid R groups and amylopectin repeat counts remain unresolved.
    def replace(rid,stoich,reason):
        r=m.reactions.get_by_id(rid);old={x.id:v for x,v in r.metabolites.items()}
        r.subtract_metabolites(dict(r.metabolites));r.add_metabolites({m.metabolites.get_by_id(k):v for k,v in stoich.items()})
        ledger.append({'reaction':rid,'before':json.dumps(old),'after':json.dumps(stoich),'reason':reason})
    # Correct the two discrepancies against published S1, then normalize the
    # amino-acid polymer residues to 1g per protein pseudo-unit.
    r=m.reactions.R610;shape={x.id:-v for x,v in r.metabolites.items() if v<0 and x.id!='ATP[c]'}
    shape.update({'LSER[c]':.357,'LLEU[c]':.492})
    water=m.metabolites.get_by_id('H2O[c]').formula_weight
    pm=sum(v*(m.metabolites.get_by_id(k).formula_weight-water) for k,v in shape.items())/1000
    shape={k:v/pm for k,v in shape.items()};g=37.195
    sto={k:-v for k,v in shape.items()};sto.update({'ATP[c]':-g,'ADP[c]':g,'Pi[c]':g,'H2O[c]':sum(shape.values())-g,'PROTEIN[c]':1})
    replace('R610',sto,'S1 serine/leucine restored; amino-acid shape normalized by residue mass. Cys/Trp and 37.195 ATP/g protein remain borrowed assumptions; neutral-formula hydrolysis and polymerization water explicit.')
    counts=collections.Counter(''.join(str(x.seq).upper() for x in SeqIO.parse(ROOT/'inputs/genomes/GCA_000283495.2/genomic.fna','fasta')))
    total=sum(counts[b] for b in 'ATGC');gc=(counts['G']+counts['C'])/total
    frac={'dATP[c]':(1-gc)/2,'dTTP[c]':(1-gc)/2,'dGTP[c]':gc/2,'dCTP[c]':gc/2}
    ppi=m.metabolites.get_by_id('PPi[c]').formula_weight
    mean=sum(v*(m.metabolites.get_by_id(k).formula_weight-ppi) for k,v in frac.items())
    coeff={k:1000*v/mean for k,v in frac.items()};g=4.403
    sto={k:-v for k,v in coeff.items()};sto.update({'ATP[c]':-g,'ADP[c]':g,'Pi[c]':g,'H2O[c]':-g,'PPi[c]':sum(coeff.values()),'DNA[c]':1})
    replace('R608',sto,'DSM20460 GCA_000283495.2 dsDNA composition; exactly 1g model-convention DNA residues. PPi equals incorporated nucleotides; inherited 4.403 ATP/g DNA remains assumed.')
    m.id='iME375_partial_DNA_protein_repair_candidate_v6'
    m.notes['limits']='Partial candidate only: undefined phospholipid R groups, amylopectin repeat count, other module water/energy and whole-cell fractions are NOT resolved. Not a calibrated primary comparator.'
    dest=ROOT/'models/iME375_partial_DNA_protein_repair_candidate_v6.xml';cobra.io.write_sbml_model(m,str(dest))
    pd.DataFrame(ledger).to_csv(OUT/'iME375_partial_repair_ledger.csv',index=False)
    record('results/comparator_calibration_v6/iME375_partial_repair_scope.json',{'source_DNA_GC':gc,'protein_residue_mass_before_normalization':pm,'new_protein_residue_mass':1.0,'new_DNA_residue_mass':1.0,'whole_cell_mass_resolved':False,'remaining':['PLIPID generic R groups','Amylopectin unspecified polymer repeat count','RNA/peptidoglycan energy and condensation conventions','E. coli-derived whole-cell fractions are not M. elsdenii measurements'],'not_fitted_to_growth_or_butyrate':True})
    return dest

def normalized_agora(src):
    audits=pd.read_csv(ROOT/'results/comparators/v4_biomass_formula_audit.csv');made={};ledger=[]
    for x in audits.itertuples():
        mass0=x.known_element_net_mass_g_per_model_unit
        if x.unresolved_elements!='{}' or mass0<=0:
            ledger.append({'model':x.model,'status':'excluded_from_mass_normalization','mass':mass0,'reason':'Uncancelled generic polymer X and nonphysical known-element net mass. No positive physical normalization available.'});continue
        m=load(path=src[x.model]);b=m.reactions.get_by_id(x.biomass_reaction)
        # Remove balanced GAM, normalize chemical mass, then put GAM back.
        hydro={m.metabolites.get_by_id(k+'[c]'):v for k,v in {'atp':-1,'h2o':-1,'adp':1,'pi':1,'h':1}.items()}
        # Formula-based mass of this hydrolysis vector must be negligible.
        dm=sum(-v*met.formula_weight for met,v in hydro.items())/1000
        assert abs(dm)<1e-6
        gam=b.get_coefficient('adp[c]');old=dict(b.metabolites);new={}
        pseudo={'dnarep[c]','proteinsynth[c]','rnatrans[c]','biomass[c]'}
        for met,v in old.items():new[met]=v if met.id in pseudo else (v-gam*hydro.get(met,0))/mass0+gam*hydro.get(met,0)
        b.subtract_metabolites(old);b.add_metabolites(new)
        delta=max(abs(new.get(k,0)-v) for k,v in old.items())
        m.id=m.id+'_formula_mass_normalized_candidate_v6';dest=ROOT/'models'/(x.model+'_formula_mass_normalized_candidate_v6.xml')
        m.notes['v6_scope']='Formula-mass normalization candidate only. Chemical biomass excluding bookkeeping pseudo-metabolites normalized to 1g; inherited GAM preserved in mmol/g. Not measured organism-specific biomass or empirical growth calibration.'
        cobra.io.write_sbml_model(m,str(dest));made[x.model+'_mass_candidate']=dest
        ledger.append({'model':x.model,'status':'formula_normalized_candidate','original_mass_g_per_unit':mass0,'new_mass_g_per_unit':1.0,'GAM_preserved':gam,'max_coefficient_change':delta,'experimental_calibration':False})
    pd.DataFrame(ledger).to_csv(OUT/'formula_mass_normalization_decisions.csv',index=False)
    return made

def cycle_candidates(partial):
    made={};ledger=[]
    for label,source in [('iME375_cycle_constrained',ROOT/'inputs/comparators/iME375/Supplementary data 4.xml'),('iME375_partial_plus_cycle',partial)]:
        m=load(path=source);r=m.reactions.R421;old=r.bounds
        assert r.get_coefficient('ATP[c]')==-2 and r.get_coefficient('ADP[c]')==2
        r.lower_bound=0;m.id=label+'_candidate_v6'
        r.notes['v6_curation']='Forward-only full carbamoyl-phosphate-synthetase reaction; reverse full reaction creates an ATP-generating cycle with carbamate kinase R331. A partial reverse ATP-synthesis activity is not the reverse of the complete two-ATP reaction.'
        r.notes['evidence']='IUBMB EC 6.3.4.16 https://iubmb.qmul.ac.uk/enzyme/EC6/3/4/16.html ; experimental CPS mechanism https://pmc.ncbi.nlm.nih.gov/articles/PMC2253338/'
        dest=ROOT/'models'/(label+'_candidate_v6.xml');cobra.io.write_sbml_model(m,str(dest));made[label]=dest
        ledger.append({'candidate':label,'reaction':'R421','before_bounds':str(old),'after_bounds':str(r.bounds),'no_other_direction_change':True,'source_sha256':sha(source),'candidate_sha256':sha(dest),'evidence':'https://iubmb.qmul.ac.uk/enzyme/EC6/3/4/16.html','reason':'Remove complete-reaction reverse CPS / carbamate-kinase ATP-generating cycle; preserve deposited source.'})
    pd.DataFrame(ledger).to_csv(OUT/'iME375_cycle_curation_ledger.csv',index=False)
    return made

def common_medium(src):
    # Build a canonical chemical inventory from deposited AGORA exchanges.
    canon={};mappings={};inventory=[]
    for name,p in src.items():
        if 'AGORA2' not in name or name.endswith('_mass_candidate'):continue
        m=load(path=p);mp={}
        for ex in m.exchanges:
            met=next(iter(ex.metabolites));key=met.id.split('[')[0]
            if key in ['biomass','o2']:continue
            seed=met.annotation.get('seed.compound','');seed=seed[0] if isinstance(seed,list) and len(seed)==1 else seed
            info={'canonical':key,'name':met.name,'carbon_atoms':met.elements.get('C',0),'seed':AA_SEED.get(key,seed),'formula':met.formula}
            if met.formula and key not in canon:canon[key]=info
            mp[key]=ex.id
        mappings[name]=mp
    # Common rich support consists of all AA and the union of previously
    # diagnosed organic growth-support compounds; sugars used as escape routes
    # are excluded. The same carbon budget includes ALL organic nutrients.
    diag=pd.read_csv(ROOT/'results/comparators/v4_minimal_extra_imports_diagnostic.csv')
    keys=set(AA_VMH)
    for row in diag.dropna(subset=['exchange']).itertuples():
        key=str(row.exchange).removeprefix('EX_').replace('(e)','').replace('[e]','')
        if key not in ['arab_L','glc_D','lac_L','lac_D']:keys.add(key)
    keys.update(AGORA.values());keys.discard('o2')
    # MJR mapping uses exact ModelSEED compounds, manual stereo-specific AA
    # canonical corrections override several legacy AGORA annotations.
    seed_to_key={str(v['seed']):k for k,v in canon.items() if isinstance(v['seed'],str) and str(v['seed']).startswith('cpd')}
    seed_to_key.update({c:k for k,c in AA_SEED.items()});seed_to_key.update(AGORA)
    for name,p in src.items():
        if name in mappings:continue
        if name.endswith('_mass_candidate'):
            mappings[name]=dict(mappings[name.removesuffix('_mass_candidate')]);continue
        m=load(path=p)
        if name.startswith('iME375'):
            mp={AGORA[c]:r for c,r in IME.items() if c in AGORA and r in m.reactions};mp.update(IME_AA)
            # iME375 has no L-lactate exchange; no stereochemical substitution.
            mappings[name]=mp
        else:
            mappings[name]={seed_to_key[met.id.split('_')[0]]:ex.id for ex in m.exchanges for met in ex.metabolites if met.id.split('_')[0] in seed_to_key}
    for key in sorted(keys):
        if key not in canon:continue
        inf=canon[key];C=inf['carbon_atoms']
        if key in ['glc_D','lac_L','lac_D','ac']:cap=0
        elif key in AA_VMH:cap=1
        elif C==0 or key=='co2':cap=1000
        elif key.startswith(('ala','gly','cys')) and key not in AA_VMH:cap=.2
        elif key in ['ocdca','ttdca','hdca','ocdcea','chsterol']:cap=.1
        else:cap=.01
        inventory.append({**inf,'uptake_capacity':cap,'source':'Explicit new common enriched comparison panel, not experimental gut diet; no model gapfilling'})
    pd.DataFrame(inventory).to_csv(OUT/'common_enriched_medium.csv',index=False)
    return inventory,mappings

def run():
    OUT.mkdir(exist_ok=True);src=files();src['iME375_partial_candidate']=ime_candidate();src.update(cycle_candidates(src['iME375_partial_candidate']));src.update(normalized_agora(src))
    inventory,mapping=common_medium(src);lookup={x['canonical']:x for x in inventory};rows=[];ledger=[];fvas=[];energy=[]
    for name,path in src.items():
        m=load(path=path);mp=mapping[name]
        if name.startswith('iME375'):m.objective=m.reactions.R617
        obj=next(iter(cobra.util.solver.linear_reaction_coefficients(m)))
        for k,info in lookup.items():ledger.append({'model':name,'canonical':k,'compound':info['name'],'mapped_exchange':mp.get(k,''),'status':'mapped' if k in mp else 'missing_no_substitute'})
        for condition in ['closed_boundaries','inorganics_only']:
            with m:
                for r in m.boundary:r.bounds=(0,0)
                forced={r.id:r.lower_bound for r in m.reactions if r.lower_bound>0}
                # ATP-cycle probes must not be infeasible merely because the
                # deposited NGAM drain still requires nutrient-derived ATP.
                for r in m.reactions:
                    if r.lower_bound>0:r.lower_bound=0
                if condition=='inorganics_only':
                    m.medium={mp[k]:100 for k in ['h2o','h','pi','co2','nh4','so4','na1','k','mg2'] if k in mp}
                if name.startswith('iME375'):
                    ids={'ATP[c]':-1,'H2O[c]':-1,'ADP[c]':1,'Pi[c]':1}
                elif name=='MJR8396C':
                    ids={k+'_c0':v for k,v in {'cpd00002':-1,'cpd00001':-1,'cpd00008':1,'cpd00009':1,'cpd00067':1}.items()}
                else:ids={k+'[c]':v for k,v in {'atp':-1,'h2o':-1,'adp':1,'pi':1,'h':1}.items()}
                atp=cobra.Reaction('V6_DIAGNOSTIC_ATP_DRAIN');atp.bounds=(0,1000)
                atp.add_metabolites({m.metabolites.get_by_id(k):v for k,v in ids.items()});m.add_reactions([atp]);m.objective=atp
                energy.append({'model':name,'condition':condition,'max_ATP_drain':m.slim_optimize(error_value=np.nan),'status':m.solver.status,'forced_lower_bounds_relaxed_for_probe_only':json.dumps(forced),'drain_mass_charge_imbalance':json.dumps(atp.check_mass_balance())})
        for substrate in ['none','glc_D','lac_L','lac_D']:
            for scale in [.5,1,2]:
                with m:
                    if substrate!='none' and substrate not in mp:
                        rows.append({'model':name,'substrate':substrate,'support_scale':scale,'status':'unmapped_stereospecific_substrate'});continue
                    medium_map={mp[k]:x['uptake_capacity']*(scale if x['carbon_atoms']>0 and k!='co2' else 1) for k,x in lookup.items() if k in mp and x['uptake_capacity']>0}
                    if substrate!='none':medium_map[mp[substrate]]=60/lookup[substrate]['carbon_atoms']
                    m.medium=medium_map
                    # Exact gross organic carbon uptake budget; secretion cannot
                    # offset imports. CO2 is inorganic and explicitly exempted.
                    coeff={}
                    for k,rid in mp.items():
                        if k not in lookup or k=='co2' or lookup[k]['carbon_atoms']==0:continue
                        ex=m.reactions.get_by_id(rid);met,v=next(iter(ex.metabolites.items()))
                        var=ex.reverse_variable if v<0 else ex.forward_variable
                        coeff[var]=lookup[k]['carbon_atoms']*abs(v)
                    con=m.problem.Constraint(0,ub=60,name='V6_common_gross_organic_C_budget');m.add_cons_vars(con);m.solver.update();con.set_linear_coefficients(coeff)
                    s=m.optimize();row={'model':name,'substrate':substrate,'support_scale':scale,'organic_C_budget':60,'status':s.status,'growth_flux_raw_model_units':s.objective_value}
                    if s.status=='optimal':
                        p=pfba(m,fraction_of_optimum=1-1e-9)
                        row['actual_gross_organic_C_uptake']=sum(v*p.raw[var.name] for var,v in coeff.items()) if hasattr(p,'raw') else sum(v*var.primal for var,v in coeff.items())
                        for k in ['glc_D','lac_L','lac_D','ac','but','ppa','h2']:
                            row[k]=float(p.fluxes.get(mp.get(k,''),np.nan))
                        if scale==1 and s.objective_value>1e-8:
                            ids=[mp[k] for k in ['lac_L','lac_D','ac','but','ppa'] if k in mp]
                            f=flux_variability_analysis(m,reaction_list=ids,fraction_of_optimum=.95,processes=1).reset_index(names='reaction')
                            f['model']=name;f['substrate']=substrate;fvas.append(f)
                    rows.append(row)
        pd.DataFrame(rows).to_csv(OUT/'common_carbon_budget_simulations.csv',index=False)
        pd.DataFrame(ledger).to_csv(OUT/'common_medium_mapping.csv',index=False)
        pd.DataFrame(energy).to_csv(OUT/'comparator_ATP_cycle_checks.csv',index=False)
        if fvas:pd.concat(fvas).to_csv(OUT/'common_exchange_FVA95.csv',index=False)
        print(name,'common carbon budget comparisons completed',flush=True)
    record('results/comparator_calibration_v6/scope.json',{'comparators_gapfilled':False,'selective_model_specific_feeding':False,'shared_organic_C_uptake_budget':60,
           'medium':'Declared equal per-compound support where represented; missing transporters retained. Support-scale sensitivity .5/1/2.',
           'interpretation':'A controlled mathematical comparator screen, not empirically calibrated physiological growth. Native, formula-normalized and partially repaired candidates reported separately.',
           'target_measured_biomass_available':False,'absolute_growth_ranking_empirically_validated':False})

if __name__=='__main__':run()
