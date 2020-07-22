# %% 
# For numerical calculations
import numpy as np
import pandas as pd
import scipy as sp
import math
import git
from scipy.integrate import odeint
from numpy import arange
from scipy.integrate import odeint
import scipy.optimize 
from scipy.optimize import leastsq
from math import exp
from collections import OrderedDict
from sklearn.linear_model import LinearRegression
pd.options.mode.chained_assignment = None

# Find home directory for repo
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir
# %% 
#Load data of the experiment where further DMSP/enzyme were added to check for enzyme degradation
# over the course of the exp.
df_add = pd.read_csv(f'{homedir}/data/raw/enz_deg/Alma1_add_exps.csv')
df_add.head()

# %% 
#We will make a fit to the data using the least squares method, to see if the 
#degradation of DMSP by Alma1 follows Michaelis-Menten kinetics

#First, we'll define a function that computes the expected concentration of DMSP
#over time if the enzyme followed Michaelis-Menten kinetics
def substrate_kinetics(so, vmax, km, time):
    '''
    Function that computes the substrate concentration over time by
    numerically integrating the recursive equation
    Parameters
    ----------
    so : float.
        Initial concentration of substrate
    vmax : float.
        Max speed of enzyme
    km : float.
        Michaelis-Menten constant of enzyme
    time : array-like.
        Time points where to evaluate function
    '''
    # Compute ∆t
    delta_t = np.diff(time)[0]
    
    # Initialize array to save substrate concentration
    substrate = np.zeros(len(time))
    
    # Modify first entry
    substrate[0] = so
    
    # Loop through time points
    for i in range(len(time[1:])):
        substrate[i+1] = substrate[i] -\
                         vmax * substrate[i] / (km + substrate[i]) * delta_t
        
    return substrate

#We will now infer V_max from the data using the substrate kinetic function:
#Define a function that computes the residuals to fit into scipy's least_squares.
def resid(vmax, so, km, time, time_exp, s_exp):
    '''
    Function that computes the residuals of the substrate concentration
    according to the numerical integration of the dynamics.
    Parameters
    ----------
    vmax : float.
        Max speed of enzyme
    so : float.
        Initial concentration of substrate
    km : float.
        Michaelis-Menten constant of enzyme
    time : array-like.
        Time points where to evaluate function
    time_exp : array-like.
        Time points where data was taken.
    s_exp : array-like.
        Experimental determination of substrate concentration
        
    Returns
    -------
    residuals of experimental and theoretical values
    '''
    # Integrate substrate concentration
    substrate = substrate_kinetics(so, vmax, km, time)
    
    # Extract substrate at experimental time points
    time_idx = np.isin(time, time_exp)
    s_theory = substrate[time_idx]
    
    return s_theory - s_exp
# %%
#Let's determine the initial V_max of the rxns, assuming that they follow Michaelis-Menten kinetics

# Filter data by experiment A (started with 5 replicates with 100 uM DMSP 
# and 1.5X Alma1, where further DMSP was added after 38 min of the start of the exp.)

df_exp_a = df_add[df_add['Experiment']=='A']
# Filter data by times less than 40 min
# This is to exclude the values after the addition of extra DMSP
df_exp_a_add_i = df_exp_a[df_exp_a['Type']=='Before']

#Group data by treatment
df_group1 = df_exp_a_add_i.groupby(['Treatment'])

# Define column names
names = ['enzyme_ul_ml_rxn', 'vmax']

# Initialize empty dataframe to save fit results
df_fit_paramls_add = pd.DataFrame(columns=names)

# Loop through enzyme concentrations
for i, (group, data) in enumerate (df_group1):

    # Define time array
    time = np.linspace(0, data.Time_min.max(), 1000)
    

    # Append experimental time points
    time_exp = data.Time_min
    time = np.sort(
        np.unique(
            np.append(time, time_exp)
        )
    )
    # Extract initial concentration
    so = data.DMSP_uM.max()
    # Extract experimental concentrations
    s_exp =  data.DMSP_uM.values
    # Define km
    km = 9000
    #Fit Vmax
    popt, _ = scipy.optimize.leastsq(
    func=resid,
    x0=100,
    args=(so, km, time, time_exp, s_exp)
    )
    vmax = popt[0]
    # Create a substrate list 
    substrate = substrate_kinetics(so, vmax, km, time)
    
    # Store parameters and group as list
    fit = (group, popt[0])  
    # Convert list to pandas Series
    series = pd.Series(fit, index=names)   
    # Append fit to dataframe
    df_fit_paramls_add = df_fit_paramls_add.append(series, ignore_index=True)

df_fit_paramls_add
# %%
#Export to csv
df_fit_paramls_add.to_csv(f'{homedir}/data/processed/enz_deg/addexp_mmfit_before_DMSP_add.csv')
# %%
#Now, we will calculate the maximum velocity after the addition of further DMSP. 
#Utilize the function to get the residuals for Alma1

# Filter data by times more than 40 min
# This is to exclude the values after the addition of extra DMSP
df_exp_a_add_f = df_exp_a[df_exp_a['Type']=='After']

#Group data by treatment
df_group2 = df_exp_a_add_f.groupby(['Treatment'])

# Define column names
names = ['enzyme_ul_ml_rxn', 'vmax']

# Initialize empty dataframe to save fit results
df_fit_paramls_add2 = pd.DataFrame(columns=names)

# Loop through enzyme concentrations
for i, (group, data) in enumerate (df_group2):

    # Define time array
    time = np.linspace(data.Time_min.min(), data.Time_min.max(), 1000)
    

    # Append experimental time points
    time_exp = data.Time_min
    time = np.sort(
        np.unique(
            np.append(time, time_exp)
        )
    )
    # Extract initial concentration
    so = data.DMSP_uM.max()
    # Extract experimental concentrations
    s_exp =  data.DMSP_uM.values
    # Define km
    km = 9000
    #Fit Vmax
    popt, _ = scipy.optimize.leastsq(
    func=resid,
    x0=100,
    args=(so, km, time, time_exp, s_exp)
    )
    vmax = popt[0]
    # Create a substrate list 
    substrate = substrate_kinetics(so, vmax, km, time)
    
    # Store parameters and group as list
    fit = (group, popt[0])  
    # Convert list to pandas Series
    series = pd.Series(fit, index=names)   
    # Append fit to dataframe
    df_fit_paramls_add2 = df_fit_paramls_add2.append(series, ignore_index=True)

df_fit_paramls_add2

# %%
#Export to csv
df_fit_paramls_add2.to_csv(f'{homedir}/data/processed/enz_deg/addexp_mmfit_after_DMSP_add.csv')
# %%
