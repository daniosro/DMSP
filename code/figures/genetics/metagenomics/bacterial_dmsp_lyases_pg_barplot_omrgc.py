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
# Load data
df_stab6_cursonetal_2018 = pd.read_csv(f'{homedir}/data/processed/genetics/stab6_cursonetal_2018_tidy.csv')
df_stab6_cursonetal_2018.head()

#Filter by percentage of genome equivalents
df_pge = df_stab6_cursonetal_2018[(df_stab6_cursonetal_2018.Variable == 'percentage_of_bacteria')]

# Remove genes that we don't want to be plotted. recA is a normalizing gene present in
# all bacteria and dsyb is an enzyme involved in DMSP synthesis.
df_pge = df_pge[(df_pge.Gene != 'RecA') 
                & (df_pge.Gene != 'DSYB') 
                & (df_pge.Gene != 'DsyB')]

# %%
#Make figure
fig = plt.figure(figsize=(2.95, 1.95), dpi=192)
#Add barplot
ax = sns.barplot(y='Number', x='Gene', data=df_pge)
#Set labels
ax.set(xlabel='Gene', ylabel='Percentage of genomes with each gene')
#Rotate x tick marks
plt.xticks(rotation=30)

#Save figure
fig.savefig(f'{homedir}/figures/genetics/metagenomics/stab6_curson2018_bact.pdf', bbox_inches='tight')

# %%
