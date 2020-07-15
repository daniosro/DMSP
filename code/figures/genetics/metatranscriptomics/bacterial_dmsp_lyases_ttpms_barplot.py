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

# %%
# Set plot style
sns.set_style("ticks")
sns.set_palette("colorblind", color_codes=True)
sns.set_context("paper")

# load data
df_enzymes = pd.read_csv(f'{homedir}/data/processed/genetics/stab8_cursonetal_2018_tidy.csv')
df_enzymes.head()

# %%
# Create dataframe for the mean transcripts per million sequences (TTPMS) for each enzyme

#Calculate the mean and standard deviation of TTPMS for each enzyme
df_averages_ttpms = df_enzymes.groupby('Enzyme').Total_transcripts_per_million_sequences.agg(['mean','std']).reset_index()
df_averages_ttpms.head()
# %%
#Make figure
fig = plt.figure(figsize=(2.95, 1.95), dpi=192)
#Add barplot
ax = sns.barplot(y='mean', x='Enzyme', data=df_averages_ttpms)
#Set labels
ax.set(xlabel='Gene', ylabel='Total transcripts per million sequences')
#Rotate x tick marks
plt.xticks(rotation=30)

#Save figure
fig.savefig(f'{homedir}/figures/genetics/metatranscriptomics/stab8_curson2018_bact.pdf', bbox_inches='tight')

# %%
