# %%
import numpy as np
import pandas as pd
import scipy as sp
import math
import matplotlib.animation as animation
from scipy.integrate import odeint
from numpy import arange
from scipy.integrate import odeint
import scipy.optimize 
from scipy.optimize import leastsq
from math import exp
from collections import OrderedDict
from sklearn.linear_model import LinearRegression
pd.options.mode.chained_assignment = None
import git

# Find home directory for repo
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir

# %% 
#Load the HPLC processed data
df_hplc= pd.read_csv(f'{homedir}/data/processed/HPLC/hplc_master_table_raw_100.csv')
# Import tables with the EA processed data
df_ea_1 = pd.read_csv(f'{homedir}/data/processed/EA/20190917_EA.csv')
# %% 
# Create empty master dataframe
df_master=pd.DataFrame()

# Fill master table
df_master['Name'] = df_ea_1['Identifier']
df_master['d34S'] = df_ea_1['Correction_of_d34S_by_true_value']
df_master['Replicate'] = df_ea_1['Replicate']
df_master['Time_min'] = (df_ea_1['Time_min']).astype(int)

#Create new ID column with the combination of enzyme, time and date
df_master['ID'] = df_master['Name'] + '_' +\
df_master['Time_min'].map(str) + '_' + df_master['Replicate'].map(str)

# Sort values
df_master = df_master.sort_values(['Name', 'Replicate', 'Time_min'])
df_master.head()
# %% 
# Append f_r (fraction of DMSP remaining) column to the master table               
df_master = df_master.merge(df_hplc.filter(['ID', 'f_R']), how='left',
                 left_on='ID', right_on='ID') 

# Append corrected DMSP concentration column to the master table  
df_master = df_master.merge(df_hplc.filter(['ID', 'Real_conc']), \
                            how='left', left_on='ID', right_on='ID')  

# Append -ln(f_r) column to the master table
df_master['minus_ln_f_R'] = abs(-np.log(df_master['f_R']))

#Append 1000 × ln (1+δ34S/1000) -> approximation to get slope
# column to the master table
df_master['d34S_approx'] = 1000 * \
np.log(1+(df_master['d34S']/1000))

df_master.head()
# %% 
# Correction of the measured d34S values, which are a combination of 
# those of DMSP and those of the cell lysate
# Filter by Alma1
df_alma1= df_master[(df_master.Name == 'Alma1')]
df_alma1.head()
# %%
# Export master table. No mixing correction was performed because the amount 
# of cell lysate in the sample was negligible.
df_alma1.to_csv(f'{homedir}/data/processed/enzymes/alma1_master.csv')
# %%
