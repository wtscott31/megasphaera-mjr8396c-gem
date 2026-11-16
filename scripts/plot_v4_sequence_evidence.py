"""Publication-resolution sequence evidence; original Figures 3-7 stay intact."""
from common import *
import os
os.environ['MPLCONFIGDIR']=str(ROOT/'logs/matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

def main():
    j=json.loads((ROOT/'results/JGI_mapping/summary.json').read_text())
    d=pd.DataFrame([{'assembly':'Original JGI','length_bp':j['original_total_bp'],'GC_percent':j['original_GC_percent_canonical'],'contigs':80,'proteins':2189},{'assembly':'Public GenBank','length_bp':j['public_total_bp'],'GC_percent':51.028979,'contigs':70,'proteins':2422}])
    d.to_csv(ROOT/'tables/S61_original_and_public_genome_statistics.csv',index=False)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,3,figsize=(11,3.7),layout='constrained');colors=['#287d8e','#8b5fbf']
    a=axes[0];a.barh(d.assembly,d.length_bp/1e6,color=colors);a.set_xlim(0,2.9);a.set_xlabel('Assembly length (Mb)');a.set_title('A  Exact segments, different assemblies',loc='left',fontsize=10)
    for i,v in enumerate(d.length_bp):a.text(v/1e6+.035,i,f'{v:,} bp',va='center',fontsize=8)
    a=axes[1];a.scatter(d.GC_percent,[0,1],s=70,c=colors);a.set_yticks([0,1],d.assembly);a.set_xlim(50.8,51.25);a.set_xlabel('GC among canonical bases (%)');a.set_title('B  Direct sequence measurements',loc='left',fontsize=10)
    for i,v in enumerate(d.GC_percent):a.annotate(f'{v:.5f}%',(v,i),xytext=(0,12),textcoords='offset points',ha='center',fontsize=9)
    a.set_ylim(-.5,1.6)
    a=axes[2];a.barh(['477 model genes'],[342],color=colors[0],label='Exact protein');a.barh(['477 model genes'],[135],left=[342],color=colors[1],label='Exact genomic ORF; different start');a.set_xlim(0,477);a.set_xlabel('Uniquely mapped genes');a.set_title('C  Public-locus mapping complete',loc='left',fontsize=10);a.text(171,0,'342',ha='center',va='center',color='white');a.text(409.5,0,'135',ha='center',va='center',color='white');a.legend(loc='lower center',bbox_to_anchor=(.5,-.64),frameon=False,fontsize=8)
    for ext in ['png','svg']:fig.savefig(ROOT/f'figures/v4_original_genome_and_gene_mapping.{ext}',dpi=600,bbox_inches='tight')
    plt.close(fig)
if __name__=='__main__':main()
