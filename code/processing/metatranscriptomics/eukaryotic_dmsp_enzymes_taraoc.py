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
df_stab3_vorobevetal_2020 = pd.read_csv(f'{homedir}/data/raw/genetics/stab3_vorobevetal_2020.csv')
df_stab3_vorobevetal_2020.head()

# %%
#Create ocean and sea regions column

#Initialize empty list
stations=[]

#Loop through rows in dataframe
for index, row in df_stab3_vorobevetal_2020.iterrows():
    # Extract Station
    Station = row.Station
    # Extract number of transcripts from total Dataframe
    if 1<= Station <=30 or Station ==153:
        stations.append('Mediterranean Sea')
    elif 128<= Station <=141 or 103<= Station <=109:
        stations.append('North Pacific Ocean')
    elif 91<= Station <=102 or 110<= Station <=127:
        stations.append('South Pacific Ocean')
    elif 80<= Station <=90:
        stations.append('Southern Ocean (S of 40S)')
    elif 66<= Station <=79:
        stations.append('South Atlantic Ocean')
    elif 31<= Station <=65:
        stations.append('Indian Ocean')
    elif 142<= Station <=155:
        stations.append('North Atlantic Ocean')
    elif 158<= Station <=208:
        stations.append('Arctic Ocean')
# Create column of Ocean and sea regions
df_stab3_vorobevetal_2020['Ocean_and_sea_regions'] = stations
df_stab3_vorobevetal_2020.head()

# %%
#Calculate reads per kilobase per million mapped reads (RPKM)
df_stab3_vorobevetal_2020['Reads_per_Kb_per_million_mapped_reads'] = df_stab3_vorobevetal_2020['Transcripts'] * 1E9
df_stab3_vorobevetal_2020.head()

# %%
# Export data table
df_stab3_vorobevetal_2020.to_csv(f'{homedir}/data/processed/genetics/df_stab3_vorobevetal_2020.csv')
# %%
