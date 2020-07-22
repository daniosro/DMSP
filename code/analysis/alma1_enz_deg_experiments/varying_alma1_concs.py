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
#Load data of DMSP degradation by different Alma1 concs.
# Load data
df = pd.read_csv(f'{homedir}/data/raw/enz_deg/Alma1_enz_deg_DMSP_100uM.csv')

# Create real concentration column
df['dmsp_um_real']= df ['dmsp_um'] * 10

#Sort values
df = df.sort_values(['enzyme_ul_ml_rxn', 'time_min'])
df.head()
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
    for i,t in enumerate(time[1:]):
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
#Let's utilize the previous function to calculate the V_max for each concentration of enzyme.
#Group data by enzyme concentration
df_group = df.groupby(['enzyme_ul_ml_rxn'])

# Define column names
names = ['enzyme_ul_ml_rxn', 'vmax']

# Initialize empty dataframe to save fit results
df_fit_paramls = pd.DataFrame(columns=names)

# Loop through enzyme concentrations
for i, (group, data) in enumerate (df_group):

    # Define time array
    time = np.linspace(0, data.time_min.max(), 1000)
    

    # Append experimental time points
    time_exp = data.time_min
    time = np.sort(
        np.unique(
            np.append(time, time_exp)
        )
    )
    # Extract initial concentration
    so = data.dmsp_um_real.max()
    # Extract experimental concentrations
    s_exp =  data.dmsp_um_real.values
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
    df_fit_paramls = df_fit_paramls.append(series, ignore_index=True)

df_fit_paramls
# %% 
#Export to csv
df_fit_paramls.to_csv(f'{homedir}/data/processed/enz_deg/alma1_varying_concs_mmfit.csv')

# %%
