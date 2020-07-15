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

# %%
# Plot boxplot with seaborn
#Make figure
fig = plt.figure(figsize=(2.95, 1.95), dpi=192)
#Make boxplot
ax=sns.boxplot(y='Total_transcripts_per_million_sequences', x='Enzyme', 
                 data=df_enzymes, 
                 width=0.6,
               linewidth=1)
#Set axes labels
ax.set(xlabel='Gene', ylabel='Total transcripts per million sequences')
#Rotate x tick marks
plt.xticks(rotation=30)
 

# Add swarmplot
# bplot=sns.swarmplot(y='Total_transcripts_per_million_sequences', x='Enzyme',
#               data=df_enzymes, 
#               color='black',
#               alpha=0.2)

#Save figure
fig.savefig(f'{homedir}/figures/genetics/metatranscriptomics/stab8_curson2018_bact_boxp.pdf', bbox_inches='tight')

# %%
