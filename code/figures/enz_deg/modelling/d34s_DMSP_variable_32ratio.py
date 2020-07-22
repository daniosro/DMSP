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
#Below, we will examine how different 32Vmax would affect the 
# 34^epsilon of DddP at a enzyme degradation rate.

#Load data
df_34e_variable_32ratio =  pd.read_csv(f'{homedir}/data/modelling/34e_variable_32ratio.csv')
df_34e_variable_32ratio.head()
# %% 
#Create figure    
fig, ax = plt.subplots(figsize=(2.95, 1.95), dpi=192)
#Plot the Vmax vs. the linear approximation to the d34s of DMSP
plt.plot(df_34e_variable_32ratio['32ratio'], df_34e_variable_32ratio['34epsilon'])

# #Adjust y axis tick marks and decimal places#ax.set_yticks(np.linspace(3.9,3.98,5))
ax.set_xticks(np.linspace(0,0.5,5))
ax.set_yticks(np.linspace(3.9,3.98,5))
plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}'))
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
# #Add title and axes labels
plt.title ('Constant enzyme deg. rate')
plt.xlabel ('$^{32}V_{max}/^{32}K_{M} (min^{-1}nM\;enzyme^{-1})$')
plt.ylabel ('$^{34}\epsilon_{DddP}$ (‰)')
ax.axvline(linewidth=1, x = 0.03, color='black', linestyle='--')

# %% 
#Save figure
fig.savefig(f'{homedir}/figures/enz_deg/modelling/change_in_34epsilon_with_32ratio.pdf', 
bbox_inches='tight')