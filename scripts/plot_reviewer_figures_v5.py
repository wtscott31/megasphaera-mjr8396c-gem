#!/usr/bin/env python3
"""Regenerate two legible manuscript figures from supplied, versioned tabular outputs.

Figure 2 uses only experimentally reported comparison values explicitly stated in
manuscript, avoiding invented acetate/propionate chemostat measurements.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from common import ROOT as RELEASE, BASE
ROOT=BASE
OUT=RELEASE/'figures'/'v5'
OUT.mkdir(parents=True,exist_ok=True)
COL={'teal':'#177f83','amber':'#e5a235','rust':'#b55b42','navy':'#2e547c','purple':'#806ba3','grey':'#68747d'}
plt.rcParams.update({'font.size':10.5,'font.family':'DejaVu Sans','axes.titlesize':12,'axes.titleweight':'bold','axes.labelsize':10,'figure.dpi':130,'savefig.dpi':600,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):
    for ext in ('pdf','svg','png'):
        fig.savefig(OUT/f'{name}.{ext}',bbox_inches='tight',pad_inches=.15)
    plt.close(fig)

md=pd.read_csv(ROOT/'tables'/'Table_3_growth_and_exchange_three_media.csv')
meds=['Gut-like','Carbohydrate-enriched','Six-amino-acid-enriched']
conditions=['Gut-like','Carbon-rich','AA-rich']
x=np.arange(3)

fig,axes=plt.subplots(2,2,figsize=(11.5,10.6),gridspec_kw={'hspace':.43,'wspace':.32})
(axA,axB),(axC,axD)=axes
axA.set_xlim(0,10);axA.set_ylim(0,9.8);axA.axis('off');axA.set_title('A   Reconstruction workflow',loc='left',pad=12)
steps=[('MJR8396C sequence inputs','Original JGI / public assembly'),('gapseq reconstruction','Pathways, transport and genes'),('Anaerobic medium + gapfilling','Documented constraints and model'),('Native numerical diagnostics','Genes, biomass, cycles and FVA'),('Controlled model scenarios','Nutrients, deletions and exchange')]
for i,(s1,s2) in enumerate(steps):
    yy=8.8-i*1.72
    axA.add_patch(FancyBboxPatch((.68,yy-1.20),8.65,1.18,boxstyle='round,pad=0.09',linewidth=1.0,edgecolor='#657587',facecolor='#f3f6f7'))
    axA.text(5.02,yy-.31,s1,ha='center',va='center',weight='bold',fontsize=10.8)
    axA.text(5.02,yy-.93,s2,ha='center',va='center',fontsize=9,color='#54616c')
    if i<4: axA.annotate('',xy=(5.02,yy-1.68),xytext=(5.02,yy-1.34),arrowprops={'arrowstyle':'-|>','color':'#657587','lw':1.35})
axB.set_title('B   Organic-acid carbon partitioning',loc='left',pad=42)
C={'Butyrate':4,'Propionate':3,'Formate':1,'Succinate':4}
vals=np.array([[max(0,float(md.iloc[i][m]))*a for m,a in C.items()] for i in range(3)])
fr=vals / vals.sum(axis=1)[:,None]*100
base=np.zeros(3)
for j,(m,a) in enumerate(C.items()):
    axB.bar(x,fr[:,j],bottom=base,width=.66,label=m,color=[COL['teal'],COL['amber'],COL['rust'],COL['purple']][j])
    for i in range(3):
        if fr[i,j]>8: axB.text(i,base[i]+fr[i,j]/2,f'{fr[i,j]:.0f}%',ha='center',va='center',fontsize=10,weight='bold',color='white' if j in (0,2,3) else '#282828')
    base+=fr[:,j]
axB.set_xticks(x,conditions);axB.set_ylim(0,108);axB.set_ylabel('Fraction of secreted organic-acid carbon (%)')
axB.legend(ncol=4,loc='lower center',bbox_to_anchor=(.5,1.015),fontsize=10,frameon=False)
axC.set_title('C   Defined nutrient-capacity changes',loc='left',pad=16)
medium=pd.read_csv(ROOT/'tables'/'Table_S1_medium_definitions_three_media.csv')
changed=medium[(medium['carbohydrate_rich']!=medium['gut_like']) | (medium['amino_acid_rich']!=medium['gut_like'])].copy()
changed=changed[changed['gut_like']>0].copy()
axC.axis('off')
cell=[]
for row in changed.itertuples():
    cell.append([row.compound_name,f'{row.carbohydrate_rich/row.gut_like:.2f}',f'{row.amino_acid_rich/row.gut_like:.2f}'])
tab=axC.table(cellText=cell,colLabels=['Nutrient','Glucose\nenriched','Six-AA\nenriched'],cellLoc='left',colWidths=[.48,.26,.26],bbox=[0,.12,1,.83])
tab.auto_set_font_size(False);tab.set_fontsize(11)
for (r,c),box in tab.get_celld().items():
    box.set_edgecolor('#d0d7dd');box.set_linewidth(.6)
    if r==0:box.set_facecolor('#edf2f5');box.set_text_props(weight='bold')
axC.text(0,.02,'Dimensionless fold change from gut-like bounds',transform=axC.transAxes,fontsize=10.5,color='#555')
axD.set_title('D   Biomass and organic-acid output',loc='left',pad=12)
growth=np.array([float(z) for z in md['Growth rate (h^-1)']]);tot=vals.sum(axis=1)
w=.34
bars=axD.bar(x-w/2,growth,width=w,color=COL['navy'],label='Max. biomass flux (h$^{-1}$)')
axD.set_ylabel('Maximum biomass flux (h$^{-1}$)',color=COL['navy']);axD.set_ylim(0,1.5)
for b in bars: axD.text(b.get_x()+b.get_width()/2,b.get_height()+.035,f'{b.get_height():.2f}',ha='center',fontsize=9)
axD2=axD.twinx();bx=axD2.bar(x+w/2,tot,width=w,color=COL['teal'],label='Secreted organic-acid carbon')
axD2.set_ylabel('Secreted organic-acid carbon (mmol C gDW$^{-1}$ h$^{-1}$)',color=COL['teal']);axD2.set_ylim(0,160)
for b in bx:axD2.text(b.get_x()+b.get_width()/2,b.get_height()+3,f'{b.get_height():.1f}',ha='center',fontsize=9)
axD.set_xticks(x,conditions)
fig.suptitle('MJR8396C genome-scale model: defined-media predictions',fontsize=15.5,weight='bold',y=.993)
save(fig,'Figure_1_reviewer_corrected')

bench=pd.read_csv(ROOT/'tables'/'Table_S5_external_benchmark_simulations.csv')
ac=bench[bench['condition'].str.startswith('Glucose + acetate')].copy()
ac['acetate_cap']=ac['condition'].str.extract(r'Glucose \+ acetate (\d+)').astype(int)
ac=ac.drop_duplicates('acetate_cap').sort_values('acetate_cap')
fig,axes=plt.subplots(2,2,figsize=(11.5,10.6),gridspec_kw={'hspace':.47,'wspace':.31})
(axA,axB),(axC,axD)=axes
axA.set_title('A   Growth-scale comparison',loc='left',pad=33)
rate=[.7004949341391764,.1561292533649124,.62]
colors=[COL['navy'],COL['navy'],COL['amber']]
b=axA.bar(range(3),rate,color=colors,width=.62)
axA.set_xticks(range(3),['MJR8396C\ngut-like','MJR8396C\nlactate proxy','M. elsdenii\nlactate batch'])
axA.set_ylabel('Maximum specific-growth flux / reported rate (h$^{-1}$)');axA.set_ylim(0,.83)
for t in b:axA.text(t.get_x()+t.get_width()/2,t.get_height()+.028,f'{t.get_height():.3f}',ha='center',weight='bold',fontsize=10)
axA.text(.0,1.02,'Different strains and media; comparison is descriptive.',transform=axA.transAxes,va='bottom',fontsize=10,color='#555')
axB.set_title('B   Product-rate comparison',loc='left',pad=10)
# Only two experimentally reported values are explicitly identified and anchored to the manuscript;
# avoid inventing measured acetate/propionate rates.
compounds=['Butyrate','CO$_2$'];pred=[11.3348,25.9781];literature=[10.4,45.6]
w=.32; xx=np.arange(2)
z=axB.bar(xx-w/2,pred,width=w,color=COL['navy'],label='MJR8396C model');z2=axB.bar(xx+w/2,literature,width=w,color=COL['amber'],label='Published M. elsdenii')
axB.set_xticks(xx,compounds);axB.set_ylabel('Product rate (mmol gDW$^{-1}$ h$^{-1}$)');axB.set_ylim(0,53)
for bar in list(z)+list(z2):axB.text(bar.get_x()+bar.get_width()/2,bar.get_height()+1.5,f'{bar.get_height():.1f}',ha='center',fontsize=9.5)
axB.legend(loc='upper left',frameon=False,fontsize=9)
# Condition matching and 39.2 uptake bound specified in caption, not overlaid on bars.
axC.set_title('C   Batch lactate organic-acid partitioning',loc='left',pad=32)
fractions=np.array([[0,57,43],[44,35,21]])
colors=[COL['navy'],COL['rust'],COL['teal']]; left=np.zeros(2)
for j,acid in enumerate(['Acetate','Propionate','Butyrate']):
    axC.barh([0,1],fractions[:,j],left=left,color=colors[j],height=.58,label=acid)
    for i in range(2):
        if fractions[i,j]>12:axC.text(left[i]+fractions[i,j]/2,i,f'{fractions[i,j]}%',ha='center',va='center',color='white',weight='bold',fontsize=10)
    left+=fractions[:,j]
axC.set_yticks([0,1],['MJR8396C model','M. elsdenii batch']);axC.set_xlim(0,100)
axC.invert_yaxis();axC.set_xlabel('Acid-product molar fraction (%)');axC.legend(ncol=3,loc='lower left',bbox_to_anchor=(0,1.02),frameon=False,fontsize=9)
axC.text(.01,-.24,'Jensen-Shannon distance = 0.52',transform=axC.transAxes,color='#555',fontsize=10.5)
axD.set_title('D   Acetate-capacity response',loc='left',pad=10)
xs=ac['acetate_cap'].to_numpy();first=ac.iloc[0]
change_g=100*(ac['growth_rate'].to_numpy()/float(first['growth_rate'])-1)
change_b=100*(ac['butyrate'].to_numpy()/float(first['butyrate'])-1)
axD.plot(xs,change_g,marker='o',ms=6,color=COL['navy'],lw=2,label='Growth change')
axD.plot(xs,change_b,marker='s',ms=6,color=COL['rust'],lw=2,label='Butyrate change')
axD.set_xlabel('Acetate uptake capacity\n(mmol gDW$^{-1}$ h$^{-1}$)')
axD.set_ylabel('Change from no-acetate model (%)');axD.set_ylim(-5,max(110,max(change_b)+8));axD.legend(loc='upper left',frameon=False,fontsize=9)

fig.suptitle('External physiological benchmarking: agreements and discrepancies',fontsize=15.5,weight='bold',y=.995)
save(fig,'Figure_2_reviewer_corrected')
print('Created',len(list(OUT.iterdir())),'files in',OUT)
