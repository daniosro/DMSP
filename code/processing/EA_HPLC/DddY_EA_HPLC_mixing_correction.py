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
df_ea_1 = pd.read_csv(f'{homedir}/data/processed/EA/20190903_EA.csv')
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
# Filter by DddY 
df_dddy= df_master[(df_master.Name == 'DddY')]
df_dddy.head()
# %%
#Conversion of delta values to isotopic ratios
df_dddy['r34_approx'] = ((df_dddy['d34S_approx']/1000)+1)*0.0450045
df_dddy.head()
# %%
#Conversion of isotopic ratios to fractional abundances
df_dddy['f34_approx'] = df_dddy['r34_approx']/(1+df_dddy['r34_approx'])
df_dddy.head()
# %%
#Create dataframe with true values for cell lysate and DMSP

# List of true values
dmsp_tv = [14.402934, 14.318284, 14.264486]
cell_lysate_tv = [7.765655, 7.971817, 8.227773] 
  
# Dictionary of lists  
dict_tv = {'d34S_dmsp': dmsp_tv, 'd34S_cell_lysate': cell_lysate_tv}  

#Send dictionary to dataframe
df_tv = pd.DataFrame(dict_tv) 
    
#Conversion of delta values to isotopic ratios
df_tv['r34_dmsp'] = ((df_tv['d34S_dmsp']/1000)+1)*0.0450045
df_tv['r34_cell_lysate'] = ((df_tv['d34S_cell_lysate']/1000)+1)*0.0450045

#Conversion of isotopic ratios to fractional abundances
df_tv['f34_dmsp'] = df_tv['r34_dmsp']/(1+df_tv['r34_dmsp'])
df_tv['f34_cell_lysate'] = df_tv['r34_cell_lysate']/(1+df_tv['r34_cell_lysate'])
# %%
#Find the concentration of S in the cell lysate
# Filter by dddy at t=0
df_dddy_t_0= df_dddy[(df_dddy.Time_min == 0)]
df_dddy_t_0

#Calculation of average values of 34_F
f34_dmsp = df_tv['f34_dmsp'].mean()
f34_cell_lysate = df_tv['f34_cell_lysate'].mean()
f34_mix_0 = df_dddy_t_0['f34_approx'].mean()

#Calculation of the concentration of S in the cell lysate in micromolar

#Define empty list
conc_s_cell_lysate= []
# Loop through the concentrations of DMSP at t=0
for i in df_dddy_t_0.Real_conc:
    #Get the concentration of s at t=0 from Eq. 3
    s = i * (f34_dmsp - f34_mix_0)/(f34_mix_0-f34_cell_lysate)
    #Append to list
    conc_s_cell_lysate.append(s)
#Append list to dataframe
df_dddy_t_0['conc_s_cell_lysate']=conc_s_cell_lysate
df_dddy_t_0
# %%
#Create dataframe just for conc of S in cell lysate in each replicate
df_s_cell_lysate = df_dddy_t_0[['Replicate', 'conc_s_cell_lysate']]
df_s_cell_lysate
# %%
# Append the concentration of S in cell lysate to the main enzyme dataframe
# for each replicate                
df_dddy = pd.merge(df_dddy, df_s_cell_lysate,  how='inner', on=['Replicate'], suffixes=('', '_y'))
df_dddy.head()
# %%
#Calculation of 34_F_DMSP
# Append column to the enzyme dataframe with the app. value of ^34F_DMSP
df_dddy['F34_approx_DMSP']= ((df_dddy['conc_s_cell_lysate'] * (df_dddy['f34_approx']-f34_cell_lysate))+ 
                             (df_dddy['Real_conc']*df_dddy['f34_approx']))/df_dddy['Real_conc']
df_dddy.head()
# %%
#Convert fractional abundances to isotopic ratios
df_dddy['R34_approx_DMSP'] = df_dddy['F34_approx_DMSP']/(1-df_dddy['F34_approx_DMSP'])
#Convert isotopic ratios to delta values
df_dddy['d34S_approx_DMSP'] = ((df_dddy['R34_approx_DMSP']/0.0450045)-1)*1000
df_dddy.head()
# %%
# Shift corrections
#Correct shift in replicate a
#Filter by replicate a
df_dddy_a= df_dddy[(df_dddy.Replicate == 'a')]
#Drop replicate b in the dddy dataframe
df_dddy = (df_dddy[df_dddy['Replicate'] != 'a'])
#Add 0.8 to d34_S_DSMP only in replicate a
df_dddy_a['d34S_approx_DMSP']= df_dddy_a['d34S_approx_DMSP'] + 0.8
#Add corrected replicate a to the dddk dataframe
df_dddy = pd.merge(df_dddy, df_dddy_a, how='outer')
df_dddy
# %%
# Export master table with the data with all corrections
df_dddy.to_csv(f'{homedir}/data/processed/enzymes/dddy_master.csv')