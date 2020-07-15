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
df_stab3_vorobevetal_2020 = pd.read_csv(f'{homedir}/data/processed/genetics/df_stab3_vorobevetal_2020.csv')
df_stab3_vorobevetal_2020.head()

# Plot boxplot with seaborn classified by ocean and sea regions
# IO = Indian Ocean, NAO = North Atlantic Ocean, 
# NPO = North Pacific Ocean, SAO = South Atlantic Ocean,
# SO = Southern Ocean, SPO = South Pacific Ocean. 

# Make boxplot for dsyb

#Filter by dsyb only
df_dsyb_euk = df_stab3_vorobevetal_2020[df_stab3_vorobevetal_2020['Gene']=='dsyb']
#Sort by ocean and sea regions
df_dsyb_euk = df_dsyb_euk.sort_values(by=['Ocean_and_sea_regions'])
# Plot boxplot with seaborn
fig = plt.figure(figsize=(2.95, 1.95), dpi=192)
ax=sns.boxplot(y='Reads_per_Kb_per_million_mapped_reads', 
               x='Ocean_and_sea_regions', 
               data=df_dsyb_euk,
                 width=0.6,
               linewidth=1)
#Set axes limits and labels
ax.set(xlabel='Ocean and sea regions', ylabel='Reads per Kb per \n million mapped reads')
ax.set(yscale="log")
group_labels = ['IO', 'MS', 'NAO', 'NPO', 'SAO', 'SPO', 'SO']
ax.set_xticklabels(group_labels)

# add swarmplot
# bplot=sns.swarmplot(y='Total_transcripts_per_million_sequences', x='Enzyme',
#               data=df_enzymes, 
#               color='black',
#               alpha=0.2)

#Remove minor ticks from y axis
plt.minorticks_off()
#Save figure
fig.savefig(f'{homedir}/figures/genetics/metatranscriptomics/transcripts_dsyb_oceanregs.pdf', bbox_inches='tight')

# %%
