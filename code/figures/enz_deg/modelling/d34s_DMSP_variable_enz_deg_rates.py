# %% 
import git
import pandas as pd
import numpy as np
from scipy.integrate import odeint

# Find home directory for repo
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir

#Import plotting features
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns

# Set plot style
sns.set_style("ticks")
sns.set_context("paper")
# %% 
#Below, we will examine how different enzyme degradation rates would affect the 
# 34^epsilon of DddP at a constant V_max.

#Load data
df_34e_variable_enz_deg_rates =  pd.read_csv(f'{homedir}/data/modelling/34e_variable_enz_deg_rates.csv')
df_34e_variable_enz_deg_rates.head()
# %% 

#Create figure    
fig, ax = plt.subplots(figsize=(2.95, 1.95), dpi=192)
#Plot the enzyme degradation rates vs. the linear approximation to the d34s of DMSP
plt.plot(df_34e_variable_enz_deg_rates['Rates_of_enz_deg'], df_34e_variable_enz_deg_rates['34epsilon'])
#slopes = np.asarray(slopes)

#Adjust y axis tick marks and decimal places
ax.set_yticks(np.linspace(3.86,3.96,5))
ax.set_xlim (0,1)
plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))

#Add title and axes labels
plt.title ('Constant $V_{max}$')
plt.xlabel ('Enzyme deg. rate ($min^{-1}$)')
plt.ylabel ('$^{34}\epsilon_{DddP}$ (‰)')

#Include a line that demarks the maximum enzyme degradation rate for any of the enzymes included in this study,
#according to the curvature of the experimental data of linear approximation to the d34s of DMSP vs. time
ax.axvline(linewidth=1, x = 0.1, color='black', linestyle='--')

# %% 
#Save figure
fig.savefig(f'{homedir}/figures/enz_deg/modelling/change_in_34epsilon_variable_enz_deg_rates.pdf', bbox_inches='tight')