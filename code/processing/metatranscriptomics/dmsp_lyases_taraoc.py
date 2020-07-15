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
#Load data for the IDs of the TARA Oceans stations
df_equiv = pd.read_csv(f'{homedir}/data/raw/genetics/tara_sts_id.csv')
df_equiv.head()
# %%
# Import the dataframe which has only the bacterial enzymes of interest (DMSP lyases), 
# with total transcripts per million sequences (TTPMS) and 
# reads per kilobase per million mapped sequences (RPKM).
df_enzymes= pd.read_csv(f'{homedir}/data/processed/genetics/stab8_cursonetal_2018_tidy.csv')
df_enzymes.head()
# %%
#Load data for the eukaryotic DMSP enzymes
df_stab3_vorobevetal_2020 = pd.read_csv(f'{homedir}/data/processed/genetics/df_stab3_vorobevetal_2020.csv')
df_stab3_vorobevetal_2020.head()

# %%
#Append a column to the transcripts dataframe with the
# ID corresponding to each station
#Initialize empty lists
ids=[]

#Loop through rows in dataframe
for index, row in df_enzymes.iterrows():
    # Extract id
    ide = row.Sample_ID
    # Extract ids from equivalence Dataframe
    eq_station = df_equiv[df_equiv.id == ide].station.iloc[0]
    # Append to list
    ids.append(eq_station)

#Append to dataframe
df_enzymes['Station'] = ids
df_enzymes.head()

# %%
#Select only the columns in the Vorobev dataframe that have corresponding entries in the bacterial dataframe
df_vorobev_filtered = df_stab3_vorobevetal_2020[df_stab3_vorobevetal_2020.Station.isin(df_enzymes['Station'])]
df_vorobev_filtered.head()

# %%
#Filter by Alma1 and drop non needed columns in the dataframe
#Filter by alma1 only
df_alma1_euk_fil = df_vorobev_filtered[df_vorobev_filtered['Gene']=='alma1']
#Eliminate unnecesary columns from alma1 dataframe
df_alma1 = df_alma1_euk_fil.drop(df_alma1_euk_fil.columns[[0, 1, 2, 3,4]], axis = 1)
#Rename column
df_alma1 = df_alma1.rename(columns={"Gene": "Enzyme"})
df_alma1.head()

# %%
#Eliminate unnecesary columns from bacterial dataframe
df_bact = df_enzymes.drop(df_enzymes.columns[[0,1,3,4,5,6,7,8,9,11,12,14]], axis = 1)
df_bact.head()

# %%
# Join the bacterial and eukaryotic dataframes together
# Join dataframes
df_total_reads = pd.concat([df_alma1, df_bact])
df_total_reads.head()

# %%
# Export data table
df_total_reads.to_csv(f'{homedir}/data/processed/genetics/DMSP_lyases.csv')

# %%
