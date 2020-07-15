# %%
import os
import numpy as np
import scipy
import scipy.optimize
import pandas as pd
import git
from collections import OrderedDict

# Find home directory for repo
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir

# Import plotting features
import matplotlib.pyplot as plt
import seaborn as sns

# Set plot style
sns.set_style("ticks")
sns.set_palette("colorblind", color_codes=True)
sns.set_context("paper")
# %%
# load data
df_total_reads = pd.read_csv(f'{homedir}/data/processed/genetics/DMSP_lyases.csv')
df_total_reads.head()

# Plot boxplot with seaborn

colors = ['#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7','#000000']
# Set your custom color palette
sns.set_palette(sns.color_palette(colors))

#Create figure
fig = plt.figure(figsize=(2.95, 1.95), dpi=192)
#Make boxplot
ax=sns.boxplot(y='Reads_per_Kb_per_million_mapped_reads', x='Enzyme', 
                 data=df_total_reads, 
                 width=0.9,
               linewidth=1, 
               color='white')

#Rotate x tick marks
plt.xticks(rotation=30)
 
# add stripplot
bplot=sns.stripplot(y='Reads_per_Kb_per_million_mapped_reads', x='Enzyme',
              data=df_total_reads)
ax.set(xlabel='Gene', ylabel='Reads per Kb per \n million mapped reads')
ax.set_yscale("symlog", linthreshy=0.1, linscaley=2)
ax.set_ylim(0,1E5)

#Save figure
fig.savefig(f'{homedir}/figures/genetics/metatranscriptomics/transcripts_dmsp_lyases.pdf', bbox_inches='tight')
# %%
