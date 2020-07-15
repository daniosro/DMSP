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

# Make boxplot for alma1

#Filter by alma1 only
df_alma1_euk = df_stab3_vorobevetal_2020[df_stab3_vorobevetal_2020['Gene']=='alma1']
#Sort by ocean and sea regions
df_alma1_euk = df_alma1_euk.sort_values(by=['Ocean_and_sea_regions'])
# Plot boxplot with seaborn
fig = plt.figure(figsize=(2.95, 1.95), dpi=192)

#Make boxplot
ax=sns.boxplot(y='Reads_per_Kb_per_million_mapped_reads', 
               x='Ocean_and_sea_regions',  
               data=df_alma1_euk,
                 width=0.9,
               linewidth=1, 
               color='white')

#Rotate x tick marks
plt.xticks(rotation=30)
 
# add stripplot
bplot=sns.stripplot(y='Reads_per_Kb_per_million_mapped_reads', x='Ocean_and_sea_regions', 
              data=df_alma1_euk, color = '#E69F00')

#Set axes labels and limits
ax.set(xlabel='Ocean and sea regions', ylabel='Reads per Kb per \n million mapped reads')
ax.set_yscale("symlog", linthreshy=0.1, linscaley=2)
ax.set_ylim(0,1E5)
group_labels = ['IO', 'MS', 'NAO', 'NPO', 'SAO', 'SPO', 'SO']
ax.set_xticklabels(group_labels)

#Remove minor ticks from y axis
plt.minorticks_off()
#Save figure
fig.savefig(f'{homedir}/figures/genetics/metatranscriptomics/transcripts_alma1_oceanregs.pdf', bbox_inches='tight')

# %%
