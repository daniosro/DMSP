# %%
import numpy as np
import pandas as pd
import git

# Find home directory for repo
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir

# Import plotting features
import matplotlib.pyplot as plt
import seaborn as sns

# %%
# Set plot style
sns.set_style("ticks")
sns.set_palette("colorblind", color_codes=True)
sns.set_context("paper")

# load data
df_enzymes = pd.read_csv(f'{homedir}/data/processed/genetics/stab8_cursonetal_2018_tidy.csv')
df_enzymes.head()

# Plot boxplot with seaborn

#Define colors
colors = ['#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7','#000000']
# Set custom color palette
sns.set_palette(sns.color_palette(colors))
#Make figure
fig = plt.figure(figsize=(2.95, 1.95), dpi=192)
#Make boxplot
ax=sns.boxplot(y='Reads_per_Kb_per_million_mapped_reads', x='Enzyme', 
                 data=df_enzymes, 
                 width=0.6,
               linewidth=1)

#Rotate x tick marks
plt.xticks(rotation=30)
 
# add swarmplot
# bplot=sns.swarmplot(y='Reads_per_Kb_per_million_mapped_reads', x='Enzyme',
#               data=df_enzymes, 
#               color='black',
#               alpha=0.2)

#Set axes labels and limits
ax.set(xlabel='Gene', ylabel='Reads per Kb per \n million mapped reads')
ax.set_ylim(-5,125)

#Save figure
fig.savefig(f'{homedir}/figures/genetics/metatranscriptomics/stab8_curson2018_bact_boxp_rpkm.pdf', bbox_inches='tight')


# %%
