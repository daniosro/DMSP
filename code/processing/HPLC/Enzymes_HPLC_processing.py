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
# Load data
df_hplc = pd.read_csv(f'{homedir}/data/raw/HPLC/hplc_master_table_raw_100.csv')
 
# Sort values
# Add real concentration column to the hplc master table
#Create empty list
real_conc=[]
#loop through rows
for i, row in df_hplc.iterrows():
    # Add real concentration to empty list if it exists
    if math.isnan(row.Calc_conc):
        real_conc.append (row.Real_conc)
    # If ther real concentration does not exist, calculate it by multiplying by 10 
    #(1:10 dilution) the calculated concentration
    else:
        real_conc.append(row.Calc_conc*10)
    
df_hplc['Real_conc'] = real_conc
# Sort values
df_hplc = df_hplc.sort_values(['Name', 'Replicate', 'Time_min'])

df_hplc.head()

# %% 
# Calculate the fraction of reactant remaining for each replicate at each time point
#Create ID column with the combination of enzyme, time and replicate
df_hplc['ID'] = df_hplc['Name'] + '_' +\
df_hplc['Time_min'].astype(int).map(str) + '_' + \
df_hplc['Replicate'].map(str)

# Create new name_date column with the combination of enzyme and replicate
df_hplc['Name_Rep'] = df_hplc['Name'] + '_' +\
df_hplc['Replicate'].map(str)

# Get the max values for corrected concentration for each enzyme and replicate and 
# append it to a new column
df_hplc['Real_conc_max'] = \
df_hplc.groupby(['Name_Rep'])['Real_conc'].transform(max)

# Get the fraction of reactant remaining for each time point
df_hplc['f_R'] = abs(np.divide(df_hplc['Real_conc'], \
               df_hplc['Real_conc_max']))
df_hplc.head()
# %%
# Export data table
df_hplc.to_csv(f'{homedir}/data/processed/HPLC/hplc_master_table_raw_100.csv')
# %%