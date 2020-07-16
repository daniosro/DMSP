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
# Import table with the repeatability data
df_hplc_qc = pd.read_csv(f'{homedir}/data/raw/HPLC/hplc_qc.csv')

#Filter by date corresponding to DddP and DmdA
df_enzyme_qc = df_hplc_qc.loc[lambda x: (x['Date'] == 20190916)| (x['Date'] == 20190920)]
# Calculate the mean and standard deviation                           
df_group_qc = df_enzyme_qc.groupby(['Date']).agg(
                      {'Calc_Conc':['mean','std']})

# Create a column for relative standard deviation
df_group_qc['rstdev'] = df_group_qc['Calc_Conc', 'std'] * 100 / df_group_qc['Calc_Conc', 'mean']
df_group_qc

# %%
