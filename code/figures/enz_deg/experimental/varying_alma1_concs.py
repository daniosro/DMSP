# %% 
import numpy as np
import pandas as pd
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
df = pd.read_csv(f'{homedir}/data/raw/enz_deg/Alma1_enz_deg_DMSP_100uM.csv')
# Create real concentration column
df ['dmsp_um_real']= df ['dmsp_um'] * 10
#Sort values
df = df.sort_values(['enzyme_ul_ml_rxn', 'time_min'])
df.head()

# %% 
# Load the Michaelis-Menten fit data
df_fit_paramls= pd.read_csv(f'{homedir}/data/processed/enz_deg/alma1_varying_concs_mmfit.csv')
df_fit_paramls
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
#Let's make the figure
# Define fig and axes
fig = plt.figure(figsize=(2.95, 1.95), dpi=192)
ax = fig.add_subplot(111)
#Group data by enzyme concentration
df_group = df.groupby(['enzyme_ul_ml_rxn'])
# Define colors
colors = sns.color_palette('colorblind', n_colors=len(df_group))
# Define markers
markers = ['o', 's', 'd', '*','^']
# Loop through replicate
for i, (group, data) in enumerate(df_group):
    # Extract initial concentration
    so = data.dmsp_um_real.max()
    # Define km
    Km = 9000
    # Extract fit vmax
    vmax = df_fit_paramls[df_fit_paramls.enzyme_ul_ml_rxn == group].vmax.values
    # Define time array
    time = np.linspace(0, data.time_min.max(), 1000)
    # Append experimental time points
    time_exp = data.time_min
    time = np.sort(
        np.unique(
            np.append(time, time_exp)
        )
    )
    # Plot fit
    ax.plot(time, substrate_kinetics(so, vmax, Km, time), c=colors[i], label="")
    
    # Plot experimental data
    ax.scatter(data.time_min, data.dmsp_um_real, color=colors[i], marker=markers[i],
               label=f"{group}X")
    #ax.set_title('DddY. Vmax fitted.')
    ax.set_ylabel(r'[DMSP] ($\mu$M)')
    ax.set_xlabel(r'Time (min)')

#Set axes limits and tick marks
ax.set_xlim(-1,40)
ax.set_ylim(-5,100)
ax.set_xticks(range(0, 50, 10))
ax.set_yticks (range(0, 110, 20))
#Set legend position
ax.legend(bbox_to_anchor=(1, 0.9), title="[Alma1]")

# %%
#Save figure
fig.savefig(f'{homedir}/figures/enz_deg/experiments/Alma1_enz_deg.pdf', bbox_inches='tight')

# %%
