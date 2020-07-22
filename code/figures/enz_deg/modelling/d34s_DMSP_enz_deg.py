# %% 
import git
import pandas as pd
import numpy as np

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
#Load modelling results
df_DMSP =  pd.read_csv (f'{homedir}/data/modelling/DMSP_enz_deg.csv')
df_DMSP.head()
#Load DddP data
# Import table with the raw data for DddP
df_dddp = pd.read_csv(f'{homedir}/data/processed/enzymes/dddp_master.csv')
df_dddp.head()
# %% 
#Define slopes (fractionation factors) to apper in the plot
#Without enzyme degradation
slope_0 = -4.04
#With enzyme degradation
slope_k = -3.98

#Create figure

fig = plt.figure(figsize=(6.1, 4), dpi=192)

ax1 = plt.subplot(221)
ax2 = plt.subplot(223)
ax3 = plt.subplot(122)

#Assing plots to axes
ax1.plot(df_DMSP['-ln_fr_no_enz_deg'], df_DMSP['d34s_no_enz_deg'], color ='#0072B2')
ax2.plot(df_DMSP['-ln_fr_enz_deg'], df_DMSP['d34s_enz_deg'], color='#D55E00')
ax3.plot(df_DMSP['Time'], df_DMSP['d34s_no_enz_deg'], 
label ='$\delta^{34}{S}$ without deg', color ='#0072B2')
ax3.plot(df_DMSP['Time'], df_DMSP['d34s_enz_deg'], 
label='$\delta^{34}{S}$ with deg', color ='#D55E00')
    
#Add axes labels and legends to the first plot
ax1.set_ylabel ('$\delta ^{34}$S DMSP$_{VCDT}$ (‰)')
ax1.set_xlabel ('-ln $f_{R}$')
ax1.set_yticks(np.linspace(14,28,5))

ax1.set_title('Fractionation of S isotopes in  \n DMSP without enz. degradation')
ax1.set_xticks(range(0,6))
ax1.set_xlim(0,5)
ax1.set_yticks(range(14,35,4))
ax1.set_ylim(14,34)
#Add epsilon
ax1.text(1, 0, '$^{34}$$\epsilon$ = %s ‰' %(-slope_0), transform=ax1.transAxes,
    verticalalignment='bottom', horizontalalignment='right', fontsize=10)

#Add axes labels and legends to the second plot
ax2.set_ylabel ('$\delta ^{34}$S DMSP$_{VCDT}$ (‰)')
ax2.set_xlabel ('-ln $f_{R}$')
ax2.set_ylim (14,19)
ax2.set_title('Fractionation of S isotopes in \n DMSP with enz. degradation')
ax2.set_xticks([0,0.2,0.4,0.6,0.8,1,1.2])
ax2.set_xlim(0,1.2)
ax2.set_yticks(range(14,20,1))
ax2.set_ylim(14,19)
#Add epsilon
ax2.text(1, 0, '$^{34}$$\epsilon$ = %s ‰' %(-slope_k), transform=ax2.transAxes,
    verticalalignment='bottom', horizontalalignment='right', fontsize=10)

#Add axes labels and legends
ax3.set_xlabel ('Time (min)')
ax3.set_ylabel ('$\delta ^{34}$S DMSP$_{VCDT}$ (‰)')
ax3.set_xlim (0,60)
ax3.set_yticks(range(14,35,4))
ax3.set_ylim (14,34)
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

#Plot the DddP data on top
# Group by replicate
df_group = df_dddp.groupby(['Replicate'])
# Define colors

# Loop through replicates
for i, (group, data) in enumerate(df_group):

#Plot experimental data
    ax3.scatter(data.Time_min, data.d34S_approx_DMSP, color = '#864410')

#Show legend
ax3.legend()
#Adjust figure dimensions and save figure
plt.tight_layout(pad=0.4, w_pad=1, h_pad=1)

# %%
#Save figure
fig.savefig(f'{homedir}/figures/enz_deg/modelling/fit_d34s_enz_deg.pdf', 
bbox_inches='tight')
# %%
