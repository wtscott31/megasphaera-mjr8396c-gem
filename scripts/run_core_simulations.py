#!/usr/bin/env python3
"""Reproduce the revised manuscript-facing MJR8396C simulations.

Requirements: numpy, pandas, scipy, lxml.
The script deliberately reads the final gapseq SBML objective rather than assuming
that a reaction named "bio1" is the optimized biomass drain.
"""
from dataclasses import dataclass
from pathlib import Path
import itertools, re
import numpy as np
import pandas as pd
from lxml import etree
from scipy import sparse
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "models" / "Megasphaera_sp_MJR8396C_filled_anaerobic.xml"
MEDIA = {
    "Gut-like": ROOT / "media" / "MJR8396C_gut_like_medium.tsv",
    "Carbon-rich": ROOT / "media" / "MJR8396C_carbon_rich_medium.tsv",
    "AA-rich": ROOT / "media" / "MJR8396C_AA_rich_medium.tsv",
}
TABLES = ROOT / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

@dataclass
class ParsedSBML:
    reaction_ids: list
    reaction_names: list
    metabolite_ids: list
    S: sparse.csr_matrix
    lb: np.ndarray
    ub: np.ndarray
    objective: np.ndarray
    objective_sense: str
    gene_rules: list

def clean_id(x, prefix):
    return x[len(prefix):] if x.startswith(prefix) else x

def parse_gapseq_sbml(path):
    tree=etree.parse(str(path)); root=tree.getroot()
    ns={'sbml':root.nsmap[None],'fbc':root.nsmap['fbc']}
    sm=root.find('sbml:model',ns); fbcns='{'+ns['fbc']+'}'
    params={}
    lop=sm.find('sbml:listOfParameters',ns)
    if lop is not None:
        for p in lop: params[p.get('id')]=float(p.get('value'))
    species=sm.find('sbml:listOfSpecies',ns)
    met_xml=[s.get('id') for s in species]
    met_ids=[clean_id(x,'M_') for x in met_xml]
    met_index={x:i for i,x in enumerate(met_xml)}
    rids=[]; names=[]; lbs=[]; ubs=[]; rows=[]; cols=[]; data=[]; gene_rules=[]
    rxns=sm.find('sbml:listOfReactions',ns)
    for j,r in enumerate(rxns):
        xrid=r.get('id'); rids.append(clean_id(xrid,'R_')); names.append(r.get('name') or '')
        lbref=r.get(fbcns+'lowerFluxBound'); ubref=r.get(fbcns+'upperFluxBound')
        lbs.append(params.get(lbref,float(lbref) if lbref and re.fullmatch(r'[-+0-9.eE]+',lbref) else 0.0))
        ubs.append(params.get(ubref,float(ubref) if ubref and re.fullmatch(r'[-+0-9.eE]+',ubref) else 1000.0))
        for tag,sign in [('sbml:listOfReactants',-1),('sbml:listOfProducts',1)]:
            part=r.find(tag,ns)
            if part is not None:
                for sr in part:
                    rows.append(met_index[sr.get('species')]); cols.append(j)
                    data.append(sign*float(sr.get('stoichiometry','1')))
        gpa=r.find('fbc:geneProductAssociation',ns); refs=[]
        if gpa is not None:
            for ref in gpa.iterfind('.//fbc:geneProductRef',ns):
                refs.append(clean_id(ref.get(fbcns+'geneProduct'),'G_'))
        gene_rules.append(' or '.join(refs))
    S=sparse.coo_matrix((data,(rows,cols)),shape=(len(met_ids),len(rids))).tocsr()
    obj=np.zeros(len(rids)); sense='maximize'
    objs=sm.find('fbc:listOfObjectives',ns)
    if objs is not None:
        active=objs.get(fbcns+'activeObjective')
        for o in objs.findall('fbc:objective',ns):
            if active is None or o.get(fbcns+'id')==active:
                sense=o.get(fbcns+'type') or 'maximize'
                lof=o.find('fbc:listOfFluxObjectives',ns)
                if lof is not None:
                    for fo in lof:
                        rid=clean_id(fo.get(fbcns+'reaction'),'R_')
                        if rid in rids: obj[rids.index(rid)]=float(fo.get(fbcns+'coefficient','1'))
                break
    return ParsedSBML(rids,names,met_ids,S,np.array(lbs),np.array(ubs),obj,sense,gene_rules)

def solve_fba(pm,lb,ub,objective=None,maximize=True,A_ub=None,b_ub=None):
    c=np.array(pm.objective if objective is None else objective,float)
    res=linprog(-c if maximize else c,A_ub=A_ub,b_ub=b_ub,A_eq=pm.S,
                b_eq=np.zeros(pm.S.shape[0]),bounds=list(zip(lb,ub)),method='highs')
    res.objective_value=float(c@res.x) if res.success else np.nan
    return res

def solve_pfba(pm,lb,ub,primary_opt):
    n=len(pm.reaction_ids); obj=pm.objective
    c=np.concatenate([np.zeros(n),np.ones(n)])
    Aeq=sparse.hstack([pm.S,sparse.csr_matrix((pm.S.shape[0],n))],format='csr')
    Aeq=sparse.vstack([Aeq,sparse.hstack([sparse.csr_matrix(obj.reshape(1,-1)),sparse.csr_matrix((1,n))])])
    beq=np.concatenate([np.zeros(pm.S.shape[0]),[primary_opt]])
    I=sparse.eye(n,format='csr')
    Aub=sparse.vstack([sparse.hstack([I,-I]),sparse.hstack([-I,-I])])
    res=linprog(c,A_ub=Aub,b_ub=np.zeros(2*n),A_eq=Aeq,b_eq=beq,
                bounds=list(zip(lb,ub))+[(0,None)]*n,method='highs')
    if res.success:
        res.fluxes=res.x[:n]; res.total_flux=float(np.sum(res.x[n:]))
    return res

pm=parse_gapseq_sbml(MODEL)
exchange_indices=[i for i,r in enumerate(pm.reaction_ids) if r.startswith('EX_')]

def find_exchange_idx(cpd):
    hits=[i for i,rid in enumerate(pm.reaction_ids) if rid.startswith('EX_') and cpd in rid]
    if not hits: return None
    return sorted(hits,key=lambda i:(not pm.reaction_ids[i].endswith('_e0'),len(pm.reaction_ids[i])))[0]

def apply_medium(df):
    lb=pm.lb.copy(); ub=pm.ub.copy()
    for i in exchange_indices: lb[i]=max(0.0,lb[i])
    mapping=[]
    for row in df.itertuples(index=False):
        i=find_exchange_idx(str(row.compounds))
        if i is None:
            mapping.append({'compound':row.compounds,'name':row.name,'reaction':'','maxFlux':row.maxFlux,'mapped':False})
        else:
            lb[i]=-float(row.maxFlux); ub[i]=max(ub[i],1000.0)
            mapping.append({'compound':row.compounds,'name':row.name,'reaction':pm.reaction_ids[i],'maxFlux':row.maxFlux,'mapped':True})
    return lb,ub,pd.DataFrame(mapping)

targets={
 'D-glucose':'cpd00027','L-lactate':'cpd00159','D-lactate':'cpd00221','acetate':'cpd00029',
 'butyrate':'cpd00211','propionate':'cpd00141','formate':'cpd00047','succinate':'cpd00036',
 'ammonia':'cpd00013','hydrogen sulfide':'cpd00239','methanethiol':'cpd00324','carbon dioxide':'cpd00011'
}
rows=[]; media={}
for name,path in MEDIA.items():
    df=pd.read_csv(path,sep='\t'); media[name]=df
    lb,ub,mapping=apply_medium(df)
    fba=solve_fba(pm,lb,ub); pf=solve_pfba(pm,lb,ub,fba.objective_value)
    row={'medium':name,'growth_rate':fba.objective_value,'total_absolute_flux':pf.total_flux,
         'mapped_medium_compounds':int(mapping.mapped.sum()),'total_medium_compounds':len(mapping)}
    for label,cpd in targets.items():
        i=find_exchange_idx(cpd); row[label+'_flux']=float(pf.fluxes[i]) if i is not None else np.nan
    rows.append(row)
pd.DataFrame(rows).to_csv(TABLES/'Table_2_revised_three_medium_simulation_summary_full.csv',index=False)

# Carbon-normalized screen.
CARBON={'D-glucose':('cpd00027',6),'D-fructose':('cpd00082',6),'D-mannose':('cpd00071',6),
        'D-ribose':('cpd00105',5),'D-xylose':('cpd00154',5),'maltose':('cpd00179',12),
        'sucrose':('cpd00076',12),'pyruvate':('cpd00020',3),'acetate':('cpd00029',2),
        'D-lactate':('cpd00221',3),'L-lactate':('cpd00159',3)}
gut=media['Gut-like'].copy()
for cpd,_ in CARBON.values():
    if cpd in set(gut.compounds): gut.loc[gut.compounds.eq(cpd),'maxFlux']=0.0
lb0,ub0,_=apply_medium(gut); bg=solve_fba(pm,lb0,ub0).objective_value
crows=[]
for label,(cpd,nc) in CARBON.items():
    lb=lb0.copy(); ub=ub0.copy(); i=find_exchange_idx(cpd)
    if i is None:
        crows.append({'substrate':label,'mapped':False}); continue
    lb[i]=-60.0/nc
    fba=solve_fba(pm,lb,ub); pf=solve_pfba(pm,lb,ub,fba.objective_value)
    row={'substrate':label,'compound':cpd,'mapped':True,'carbon_atoms':nc,'uptake_bound':60/nc,
         'background_growth':bg,'growth_rate':fba.objective_value,
         'incremental_growth':fba.objective_value-bg,
         'incremental_growth_per_C':(fba.objective_value-bg)/60}
    for product,pcpd in {'butyrate':'cpd00211','propionate':'cpd00141','formate':'cpd00047','succinate':'cpd00036'}.items():
        j=find_exchange_idx(pcpd); row[product+'_flux']=pf.fluxes[j] if j is not None else np.nan
        row[product+'_yield_per_C']=row[product+'_flux']/60
    crows.append(row)
pd.DataFrame(crows).to_csv(TABLES/'Table_S2_carbon_normalized_substrate_screen.csv',index=False)

# Leave-one-out.
aa={'L-histidine':'cpd00119','L-methionine':'cpd00060','L-cysteine':'cpd00084',
    'L-aspartate':'cpd00041','glycine':'cpd00033','L-lysine':'cpd00039'}
gut_full=media['Gut-like']
base_growth=pd.DataFrame(rows).query("medium=='Gut-like'").growth_rate.iloc[0]
out=[]
for label,cpd in aa.items():
    d=gut_full.copy(); d.loc[d.compounds.eq(cpd),'maxFlux']=0.0
    lb,ub,_=apply_medium(d); r=solve_fba(pm,lb,ub)
    out.append({'removed':label,'growth_rate':r.objective_value,'relative_growth':r.objective_value/base_growth})
pd.DataFrame(out).to_csv(TABLES/'Table_S8_amino_acid_leave_one_out_revised.csv',index=False)

# Pair interactions.
aa_bg=gut_full.copy(); aa_bg.loc[aa_bg.compounds.isin(aa.values()),'maxFlux']=0.0
lb,ub,_=apply_medium(aa_bg); mu0=solve_fba(pm,lb,ub).objective_value
single={}
for label,cpd in aa.items():
    d=aa_bg.copy(); d.loc[d.compounds.eq(cpd),'maxFlux']=2.0
    lb,ub,_=apply_medium(d); single[label]=solve_fba(pm,lb,ub).objective_value
pairs=[]
for (a,ca),(b,cb) in itertools.combinations(aa.items(),2):
    d=aa_bg.copy(); d.loc[d.compounds.eq(ca),'maxFlux']=2.0; d.loc[d.compounds.eq(cb),'maxFlux']=2.0
    lb,ub,_=apply_medium(d); fba=solve_fba(pm,lb,ub); pf=solve_pfba(pm,lb,ub,fba.objective_value)
    row={'amino_acid_1':a,'amino_acid_2':b,'pair':a+' + '+b,'growth_rate':fba.objective_value,
         'interaction_effect':fba.objective_value-single[a]-single[b]+mu0}
    for product,pcpd in {'butyrate':'cpd00211','propionate':'cpd00141','formate':'cpd00047',
                          'ammonia':'cpd00013','hydrogen sulfide':'cpd00239','methanethiol':'cpd00324'}.items():
        j=find_exchange_idx(pcpd); row[product+'_flux']=pf.fluxes[j] if j is not None else np.nan
    pairs.append(row)
pd.DataFrame(pairs).to_csv(TABLES/'Table_S3_amino_acid_pair_simulations.csv',index=False)

# Butyrate envelope.
lb,ub,_=apply_medium(gut_full); max_growth=solve_fba(pm,lb,ub).objective_value
biomass_i=int(np.flatnonzero(pm.objective)[0]); but_i=find_exchange_idx('cpd00211')
env=[]
for fraction in np.linspace(0,1,31):
    l=lb.copy(); l[biomass_i]=fraction*max_growth
    obj=np.zeros(len(pm.reaction_ids)); obj[but_i]=1
    mn=solve_fba(pm,l,ub,obj,maximize=False); mx=solve_fba(pm,l,ub,obj,maximize=True)
    env.append({'growth_fraction':fraction,'growth_rate':fraction*max_growth,
                'butyrate_min':mn.objective_value,'butyrate_max':mx.objective_value})
pd.DataFrame(env).to_csv(TABLES/'Table_S4_butyrate_production_envelope.csv',index=False)
print("Finished revised simulations.")
