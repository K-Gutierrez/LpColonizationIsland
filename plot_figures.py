import config
import pandas as pd
from scipy.optimize import curve_fit
from matplotlib import cm
from scipy.stats import poisson
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.patches as mpatches

matplotlib.rc('axes', linewidth=.5)
matplotlib.rc('font', weight='bold')

label_size=20
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size 
matplotlib.rcParams['axes.spines.right'] = False
matplotlib.rcParams['axes.spines.top'] = False

color_mood = ["#a63a33","#f1c132","#3c414d","#d7d7de","#aaaaaa"]

font = {'family': 'sans-serif',
        'color':  'k',
        'size': 40,
        "weight":"bold"
        }

class fitClass:

    def __init__(self):
        pass

    def sigmoid(self, t, s):
        return self.p0*np.exp(s*t)/(1 + self.p0*(np.exp(s*t) - 1))

## timepoints
xdata = np.array([0,5,9,16,23,30,37,44,51])
rep3 = config.rep3
rep3_dic = config.rep3_dic


## read in genome features
base_genes_dir = config.base_genes_dir
df_gf = pd.read_csv(f"{config.midas_dir}/Lactobacillus/rep_genomes/Lactobacillus/genome.features",sep="\t",index_col=0)
df_gf = df_gf.loc[df_gf["gene_type"] == "CDS"]
df_gf = df_gf.loc[[d for d in df_gf.index if d in rep_genome.index]]


### Filter out genes which have ANI > 95% with any other gene 
centroids = pd.read_csv(f"{config.midas_dir}/Lactobacillus/pan_genomes/Lactobacillus/gene_info.txt",sep="\t",index_col=0)
centroids.index=centroids.index.str.replace("fig|","")
centroids = centroids.loc[df_gf.index]

for C in centroids.columns[2:]:
    centroids[C] = centroids[C].str.replace("fig|","")
    
cent = "centroid_95"

cond1 = (centroids.index != centroids[cent])
cond2 = ([c in df_gf.index for c in centroids[cent].values])
cond3 = ([c in df_gf.index for c in centroids.index])
conds=cond1*cond2*cond3

redundant_genes = centroids.loc[conds]
redundant_genes = list(set(list(redundant_genes.index) + list(redundant_genes[cent])))


## read in gene copy numbers
dfg = pd.read_csv(f"{base_genes_dir}/genes_copynum.txt",sep="\t",index_col=0)
dfg.index = [d.split("|")[1] for d in dfg.index]
dfg = dfg[rep3]
dfg = dfg.loc[[e for e in dfg.index if e in df_gf.index and e not in redundant_genes]]

## read in colonization island genes
CI_genes = pd.read_csv("colonization_island.gff",sep="\t")
CI_genes = [c.split(";")[0].split("|")[1] for c in CI_genes["##gff-version 3"].values]

## all genes on pKG
pKG_genes = dfg.loc[[f for f in df_gf.groupby("scaffold_id").get_group("contig_2").index if f in dfg.index]]
## CI genes
lost_genes = dfg.loc[[c for c in CI_genes if c in dfg.index]]
## non-pKG genes
non_lost_genes = pKG_genes.loc[[l for l in pKG_genes.index if l not in lost_genes.index]]

## contig 3 is the main chromosome
contigs = df_gf.loc[dfg.index]["scaffold_id"]
core_genome_genes = dfg.loc[[f for f in df_gf.groupby("scaffold_id").get_group("contig_3").index if f in dfg.index]]

## Plot figure showing gene copy number dynamics
fig,ax = plt.subplots(figsize=(16,8))

ax.set_facecolor("white")
ax.grid(False)

ax.plot(xdata,core_genome_genes.T,color="grey",alpha=.01,label="Other genes")

ax.plot(xdata,non_lost_genes.T.values,color="blue",alpha=.5,label="pKG genes")

ax.plot(xdata,lost_genes.T.values,color="tomato",alpha=.75)

ax.set_xlabel("Time (days)",size=25)
ax.set_ylabel("Gene copy number",size=25)

CI_patch = mpatches.Patch(color='tomato', label='Colonization island genes')
pkg_patch = mpatches.Patch(color='blue', label='pKG-WF\n(non-colonization island)')
core_patch = mpatches.Patch(color='grey', label='Core genome')

fig.legend(handles=[CI_patch,pkg_patch,core_patch],prop={"size":20},bbox_to_anchor=(.9,1.03))

ax.semilogy()

fig.savefig("figures/gene_copy_number_longitudinal", bbox_inches='tight');

## Calculate the median copy number among CI genes vs. other genes on pKG-WF
lost_genes_median = lost_genes.median()/(pKG_genes.loc[[l for l in pKG_genes.index if l not in lost_genes.index]].median())
lost_genes_median.index = pd.Series(rep3_dic).loc[lost_genes_median.index]

## estimate of 5 generations/day
gen_per_day = 5

## fit ratio of copy numbers to sigmoid function
ydata_CI = lost_genes_median.values
inst_CI = fitClass()

inst_CI.p0 = ydata_CI[0]
coeffs, coeffs_cov = curve_fit(inst_CI.sigmoid, gen_per_day*xdata, ydata_CI,-0.1)
s_CI = coeffs[0] 

t_range_gens = np.linspace(0,51*gen_per_day,10000)
sig_pred = inst_CI.sigmoid(t_range_gens,s_CI)

## plot figure
fig,ax = plt.subplots(figsize=(12,8))

ax.plot(xdata,lost_genes_median.values,lw=10,color="#2DB44A",
        marker="o",markersize=20,markeredgecolor="k",
        label="Colonization island/pKG-WF")

fit_eq = r"$f(t) = \frac{f_0e^{{s}t}}{1+f_0(e^{st} - 1)}$" 
l2 = fr"$f_0 \approx {np.around(inst_CI.p0,2)}$"
l3 = fr"$s \approx {np.around(s_CI,2)}$"
label1 = "\n" + fit_eq + "\n" + l2 + ", " + l3 

t_range = np.linspace(0,51,10000)
ax.plot(t_range,sig_pred,ls="--",lw=5,color="#131515",label=label1)

ax.set_ylabel("Frequency",size=25)
ax.set_xlabel("Time (days)",size=25)

fig.legend(prop={"size":20},bbox_to_anchor=(.95,.9))

fig.savefig("figures/CI_frequency_fit", bbox_inches='tight');
