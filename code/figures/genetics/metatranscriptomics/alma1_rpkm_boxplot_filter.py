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

# Set plot style
sns.set_style("ticks")
sns.set_palette("colorblind", color_codes=True)
sns.set_context("paper")
# %%
# load data
df_stab3_vorobevetal_2020 = pd.read_csv(f'{homedir}/data/processed/genetics/df_stab3_vorobevetal_2020.csv')
df_stab3_vorobevetal_2020.head()

# Plot boxplot with seaborn classified by filter sizes

#Filter by alma1 only
df_alma1_euk = df_stab3_vorobevetal_2020[df_stab3_vorobevetal_2020['Gene']=='alma1']

#Sort by filter
df_alma1_euk_f = df_alma1_euk.sort_values(by=['Filter'])
# Plot boxplot with seaborn
fig = plt.figure(figsize=(2.95, 1.95), dpi=192)
ax=sns.boxplot(y='Reads_per_Kb_per_million_mapped_reads', 
               x='Filter',  
               data=df_alma1_euk,
               palette=["#0072B2"],
                 width=0.6,
               linewidth=1)

# add swarmplot
# ax=sns.swarmplot(y='Reads_per_Kb_per_million_mapped_reads', x='Ocean_and_sea_regions',
#               data=df_alma1_euk, 
#               color='black',
#               alpha=0.2)

#Set axes limits and labels
ax.set(xlabel='Filter size $(\mu m)$', ylabel='Reads per Kb per \n million mapped reads')
ax.set(yscale="log")
ax.set_ylim(1,1E5)
ax.set_yticks(np.geomspace(1,1E5,6))

#Remove minor ticks from y axis
plt.minorticks_off()

#Save figure
fig.savefig(f'{homedir}/figures/genetics/metatranscriptomics/transcripts_alma1_filter.pdf', bbox_inches='tight')

# %%
