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

# Create a dataframe with the averages of RPKM for each 
# ocean and sea region
# Get averages for the reads per Kb per million mapped reads for each ocean and sea region and enzyme
df_averages_rpkm = df_enzymes.groupby(['Ocean_and_sea_regions', 'Enzyme']).Reads_per_Kb_per_million_mapped_reads.agg(['mean','std']).reset_index()
#Pivot table  to classify the means by oceanic region and enzyme
df_averages_rpkm = df_averages_rpkm.pivot_table(index=['Ocean_and_sea_regions'], columns='Enzyme', values='mean').reset_index()
#Set new index
df_averages_rpkm = df_averages_rpkm.set_index('Ocean_and_sea_regions')
#Reorder rows
df_averages_rpkm = df_averages_rpkm.reindex(['IO', 'NAO', 'NPO','SAO', 'SPO', 'SO'])
df_averages_rpkm

# %%
# Make a stacked barplot for RPKM of the bacterial DMSP degrading enzymes across ocean and sea regions.
#Let's extract the unique names of ocean and sea regions
ocean_and_sea_regions = df_enzymes[('Ocean_and_sea_regions')].unique()
# Define the x locations for the groups
N = len(ocean_and_sea_regions)
ind = np.arange(N)
#Define the width of the bars
width=0.5

#Create figure
fig = plt.figure(figsize=(2.95, 1.95), dpi=192)
ax = fig.add_subplot(111)

#Plot data
ax.bar(ind, df_averages_rpkm['alma1'], width, label='alma1', color='#000000')
ax.bar(ind, df_averages_rpkm['dddD'], width,
             bottom=df_averages_rpkm['alma1'], label='dddD', color='#56B4E9')
ax.bar(ind, df_averages_rpkm['dddK'], width,
             bottom=np.array(df_averages_rpkm['dddD'])+np.array(df_averages_rpkm['alma1']), 
             label='dddK', color = '#009E73')
ax.bar(ind, df_averages_rpkm['dddL'], width,
             bottom=np.array(df_averages_rpkm['dddK'])+np.array(df_averages_rpkm['dddD'])+
             np.array(df_averages_rpkm['alma1']), label='dddL', color='#F0E442')
ax.bar(ind, df_averages_rpkm['dddP'], width,
             bottom=np.array(df_averages_rpkm['dddL'])+np.array(df_averages_rpkm['dddK'])+
             np.array(df_averages_rpkm['dddD'])+np.array(df_averages_rpkm['alma1']), 
             label='dddP', color = '#0072B2')
ax.bar(ind, df_averages_rpkm['dddQ'], width,
             bottom=np.array(df_averages_rpkm['dddL'])+np.array(df_averages_rpkm['dddK'])+
             np.array(df_averages_rpkm['dddD'])+np.array(df_averages_rpkm['alma1'])+
            np.array(df_averages_rpkm['dddP']), 
             label='dddQ', color= '#D55E00')
ax.bar(ind, df_averages_rpkm['dddW'], width,
             bottom=np.array(df_averages_rpkm['dddL'])+np.array(df_averages_rpkm['dddK'])+
             np.array(df_averages_rpkm['dddD'])+np.array(df_averages_rpkm['alma1'])+
            np.array(df_averages_rpkm['dddP'])+np.array(df_averages_rpkm['dddQ']), 
             label='dddW',color='#D55E00')
ax.bar(ind, df_averages_rpkm['dddY'], width,
             bottom=np.array(df_averages_rpkm['dddL'])+np.array(df_averages_rpkm['dddK'])+
             np.array(df_averages_rpkm['dddD'])+np.array(df_averages_rpkm['alma1'])+
            np.array(df_averages_rpkm['dddP'])+np.array(df_averages_rpkm['dddQ'])+
                     np.array(df_averages_rpkm['dddW']), 
             label='dddY', color='#CC79A7')

#Add the labels of the x axis
plt.xticks(ind, ('IO', 'NAO', 'NPO', 'SAO', 'SPO','SO'))
ax.set(xlabel='Ocean and sea regions', ylabel='Reads per Kb per \n million mapped reads')
#Set legend
ax.legend(bbox_to_anchor=(1.2, -0.3), ncol=4)
#Set axes limits
ax.set_ylim(0,120)

#Save figure
fig.savefig(f'{homedir}/figures/genetics/metatranscriptomics/stab8_curson2018_bact_boxp_rpkm_oceanregs.pdf', bbox_inches='tight')

# %%
