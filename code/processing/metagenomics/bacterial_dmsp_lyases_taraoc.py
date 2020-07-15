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
df_stab3_landaetal_2019 = pd.read_csv(f'{homedir}/data/raw/genetics/stab3_landaetal_2019.csv')
df_stab3_landaetal_2019.head()
# %%
#Tidy up the dataframe
#Create new dataframe that will keep the first 7 columns of the old dataframe and that will have two new columns.
#Var name is the name of the categorical variable and Value name is the number assigned to each of those categories.
df_stab3_landaetal_2019 = pd.melt(
    df_stab3_landaetal_2019, 
    id_vars=df_stab3_landaetal_2019.columns.tolist()[:7],
value_vars=df_stab3_landaetal_2019.columns.tolist()[7:],
var_name='Gene', 
value_name='Percentage_of_genome_equivalents')
#Fill out the column of the categorical variable. In this case, the type of the values is not changing
df_stab3_landaetal_2019['Gene'] = df_stab3_landaetal_2019['Gene'].str[:]#.astype(str)
#Sort values
df_stab3_landaetal_2019 = df_stab3_landaetal_2019.sort_values(by=['Experiment', 'Gene'])
df_stab3_landaetal_2019.head()
# %%
# Export the dataframe
df_stab3_landaetal_2019.to_csv(f'{homedir}/data/processed/genetics/stab3_landaetal_2019_tidy.csv')

# %%
