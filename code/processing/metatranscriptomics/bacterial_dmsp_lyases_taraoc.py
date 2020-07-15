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
# %%
#Load data
df_stab8_cursonetal_2018 = pd.read_csv(f'{homedir}/data/raw/genetics/stab8_cursonetal_2018.csv')
df_stab8_cursonetal_2018.head()
# %%
# Calculate and add total transcript per million sequences (TTPMS) column
#Filter by total
df_total = df_stab8_cursonetal_2018[(df_stab8_cursonetal_2018.Enzyme == 'Total')]

#Get the total transcript number
transcript_number = df_total.Number_of_transcripts.sum()

#Initialize empty lists
transcripts=[]
ttpms = []

#Loop through rows in dataframe
for index, row in df_stab8_cursonetal_2018.iterrows():
    # Extract Sample_ID 
    Sample_ID = row.Sample_ID
    # Extract number of transcripts from total Dataframe
    transcript = df_total[df_total.Sample_ID == Sample_ID].Number_of_transcripts.iloc[0]
    # Compute Total_transcripts_per_million_sequences
    transcript_ratio = row.Number_of_transcripts * 1E6 / transcript 
    # Append to list
    ttpms.append(transcript_ratio)

# Create column of total_transcripts_per_million_sequences in dataframe
df_stab8_cursonetal_2018['Total_transcripts_per_million_sequences'] = ttpms
#Replace Alma1 for small initials
df_stab8_cursonetal_2018['Enzyme'] = df_stab8_cursonetal_2018['Enzyme'].replace(['Alma1'], 'alma1')
df_stab8_cursonetal_2018.head()

# %%
# Create dataframe for gene length in base pairs
gene_length = {'Gene':['alma1', 'dddD', 'dddK', 'dddL', 'dddP', 'dddW', 'dddY', 'dddQ', 'dmdA', 'dsyb'], 
               'Length':[1380, 2481, 756, 672, 1260, 465, 1203, 576, 1140,1050]}
df_glength = pd.DataFrame(gene_length)
df_glength

# %%
#Filter to include only enzymes
df_enzymes = df_stab8_cursonetal_2018[(df_stab8_cursonetal_2018.Enzyme != 'Total')
                                     & (df_stab8_cursonetal_2018.Enzyme != 'dsyB')
                                     & (df_stab8_cursonetal_2018.Enzyme != 'DSYB')]
df_enzymes.head()

# %%
#Calculate and add reads per kilo base per million mapped reads (RPKM) column
#Initialize empty list
rpkm_=[]

#Loop through rows in dataframe
for index, row in df_enzymes.iterrows():
    # Extract Enzyme name 
    Enzyme = row.Enzyme
    # Extract lengths from gene length Dataframe
    gl = df_glength[df_glength.Gene == Enzyme].Length.values[0]
    # Compute RPKM (Reads per kilo base per million mapped reads)
    rpkm = row.Total_transcripts_per_million_sequences * 1E3 / gl
    # Append to list
    rpkm_.append(rpkm)
# Create column of reads per kilo base per million mapped reads in dataframe
df_enzymes = df_enzymes.assign(
    Reads_per_Kb_per_million_mapped_reads = rpkm_
)

df_enzymes.head()

# %%
#Export tidy dataframe
df_enzymes.to_csv(f'{homedir}/data/processed/genetics/stab8_cursonetal_2018_tidy.csv')

# %%
