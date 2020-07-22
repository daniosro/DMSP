# %% 
import numpy as np
import pandas as pd
import scipy as sp
import scipy.optimize 
from scipy.optimize import leastsq
import git

# Find home directory for repo
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir

# Import plotting features
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.animation as animation
import seaborn as sns

# Set plot style
sns.set_style("ticks")
sns.set_palette("colorblind", color_codes=True)
sns.set_context("paper")
# %% 
#Load the experimental data
df_add = pd.read_csv(f'{homedir}/data/raw/enz_deg/Alma1_add_exps.csv')
df_add.head()

# %% 
# Load the Michaelis-Menten fit data before extra DMSP addition
df_fit_paramls_add_b= pd.read_csv(f'{homedir}/data/processed/enz_deg/addexp_mmfit_before_alma1_add.csv')
df_fit_paramls_add_b
# %%  
# Load the Michaelis-Menten fit data after extra DMSP addition
df_fit_paramls_addb2= pd.read_csv(f'{homedir}/data/processed/enz_deg/addexp_mmfit_after_alma1_add.csv')
df_fit_paramls_addb2
# %% 
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
# %% 
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
#Let's make the figure
# Define fig and axes
fig = plt.figure(figsize=(2.95, 1.95), dpi=192)
ax = fig.add_subplot(111)

# Filter data by experiment B (further addition of DMSP)
df_exp_b = df_add[df_add['Experiment']=='B']

#Group data by treatment to plot all data as scatter
df_group = df_exp_b.groupby(['Treatment'])

# Define colors
colors = sns.color_palette('colorblind', n_colors=len(df_group))
# Define markers
markers = ['o', 's', 'd', '*','^']

# Filter data by times less than 40 min
# This is to exclude the values before the addition of extra enzyme
df_exp_b_add_i = df_exp_b[df_exp_b['Type']=='Before']

# Filter data by times more than 40 min
# This is to exclude the values after the addition of extra enzyme
df_exp_b_add_f = df_exp_b[df_exp_b['Time_min']>36]

#Group data by treatment to plot all data as scatter
df_groupb = df_exp_b.groupby(['Treatment'])
#Group data before the addition of alma1 by treatment to plot the fit on top of the data
df_group_ib = df_exp_b_add_i.groupby(['Treatment'])
#Group data after the addition of alma1 by treatment to plot the fit on top of the data
df_group_fb = df_exp_b_add_f.groupby(['Treatment'])

#Generate the fit for the data before the addition of enzyme
# Loop through replicate
for i, (group, data) in enumerate(df_group_ib):
    # Extract initial concentration
    so = data.DMSP_uM.max()
    # Extract km
    Km = 9000
    # Extract fit vmax
    vmax = df_fit_paramls_add_b[df_fit_paramls_add_b.enzyme_ul_ml_rxn == group].vmax.values
    # Define time array
    time = np.linspace(0, data.Time_min.max(), 1000)
    # Append experimental time points
    time_exp = data.Time_min
    time = np.sort(
        np.unique(
            np.append(time, time_exp)
        )
    )
    # Plot fit
    ax.plot(time, substrate_kinetics(so, vmax, Km, time), c=colors[i], label="")
    
#Generate the fit for the data after the addition of enzyme
# Loop through enzyme concentrations
for i, (group, data) in enumerate (df_group_fb):

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
    
    # Plot fit
    ax.plot(time, substrate_kinetics(so, vmax, Km, time), c=colors[i], label="")

# Define labels for plots
labels = ('X','2X','3X','6X','10X')
#Loop through all data to plot them as scatter
for i, (group, data) in enumerate(df_groupb):
    # Plot experimental data
    ax.scatter(data.Time_min, data.DMSP_uM, color=colors[i], marker=markers[i],
               label=labels[i])

#Set axes labels, limits and tick marks
ax.set_ylabel(r'[DMSP] ($\mu$M)')
ax.set_xlabel(r'Time (min)')
ax.set_xlim(-1,80)
ax.set_ylim(-5,100)
ax.set_xticks(range(0, 90, 10))
ax.set_yticks (range(0, 110, 20))

# Set vertical dashed line
ax.axvline(linewidth=1, x = 38, color='black', linestyle='--')
# Set legend position
ax.legend(bbox_to_anchor=(1, -0.3), title="[Alma1]", ncol=3) 

# %%
#Save figure
fig.savefig(f'{homedir}/figures/enz_deg/experiments/Alma1_enz_deg_further_Alma1.pdf', bbox_inches='tight')
# %%
