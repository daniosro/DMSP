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
df_stab6_cursonetal_2018 = pd.read_csv(f'{homedir}/data/raw/genetics/stab6_cursonetal_2018.csv')
df_stab6_cursonetal_2018

# %%
#Tidy up the dataframe
#Create new dataframe that will keep the first 7 columns of the old dataframe and that will have two new columns.
#Var name is the name of the categorical variable and Value vars is the number assigned to each of those categories.
df_stab6_cursonetal_2018 = pd.melt(
    df_stab6_cursonetal_2018, 
    id_vars=df_stab6_cursonetal_2018.columns.tolist()[:1],
value_vars=df_stab6_cursonetal_2018.columns.tolist()[1:],
var_name='Gene')
df_stab6_cursonetal_2018.columns = ['Variable','Gene','Number']
#Sort values
df_stab6_cursonetal_2018 = df_stab6_cursonetal_2018.sort_values(by=['Variable'])

#Replace the names of the genes for small initials
df_stab6_cursonetal_2018['Gene'] = df_stab6_cursonetal_2018['Gene'].replace(['Alma1'], 'alma1')
df_stab6_cursonetal_2018['Gene'] = df_stab6_cursonetal_2018['Gene'].replace(['DddY'], 'dddy')
df_stab6_cursonetal_2018['Gene'] = df_stab6_cursonetal_2018['Gene'].replace(['DddP'], 'dddp')
df_stab6_cursonetal_2018['Gene'] = df_stab6_cursonetal_2018['Gene'].replace(['DddW'], 'dddw')
df_stab6_cursonetal_2018['Gene'] = df_stab6_cursonetal_2018['Gene'].replace(['DddK'], 'dddk')
df_stab6_cursonetal_2018['Gene'] = df_stab6_cursonetal_2018['Gene'].replace(['DddQ'], 'dddq')
df_stab6_cursonetal_2018['Gene'] = df_stab6_cursonetal_2018['Gene'].replace(['DddL'], 'dddl')
df_stab6_cursonetal_2018['Gene'] = df_stab6_cursonetal_2018['Gene'].replace(['DddD'], 'dddd')

df_stab6_cursonetal_2018.head()
 # %%
#Export the tidy dataframe
df_stab6_cursonetal_2018.to_csv(f'{homedir}/data/processed/genetics/stab6_cursonetal_2018_tidy.csv')

# %%