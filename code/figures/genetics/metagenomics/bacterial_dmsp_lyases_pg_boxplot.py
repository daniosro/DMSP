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

# Load data
df_stab3_landaetal_2019 = pd.read_csv(f'{homedir}/data/processed/genetics/stab3_landaetal_2019_tidy.csv')
df_stab3_landaetal_2019.head()

#Filter out the recA hits
df_stab3_landaetal_2019_enzymes = df_stab3_landaetal_2019[df_stab3_landaetal_2019.Gene != 'recA_hits'] 

# Plot boxplot with seaborn
fig = plt.figure(figsize=(2.95, 1.95), dpi=192)
ax=sns.boxplot(y='Percentage_of_genome_equivalents', x='Gene', 
                 data=df_stab3_landaetal_2019_enzymes, 
                 width=0.6,
                 palette="husl",
               linewidth=1)
ax.set(xlabel='Gene', ylabel='Percentage of genomes with each gene')
ax.set(ylim = (0,80))
 
# add swarmplot
# bplot=sns.swarmplot(y='Percentage_of_genome_equivalents', x='Gene',
#               data=df_stab3_landaetal_2019_enzymes, 
#               color='black',
#               alpha=0.75)
#Save figure
fig.savefig(f'{homedir}/figures/genetics/metagenomics/stab3_landa2019_bact_boxp.pdf', bbox_inches='tight')
# %%
