# %% 
import git
import pandas as pd

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
sns.set_palette("colorblind", color_codes=True)
sns.set_context("paper")

 # %% 
#Load data
df_DMSP =  pd.read_csv(f'{homedir}/data/modelling/DMSP_enz_deg.csv')
df_DMSP.head()
# %% 
#Create figure
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(6, 2), dpi=192)

#Create first plot (left; change in 34DMSP, 32DMSP and total DMSP over time)
ax1.plot(df_DMSP['Time'], df_DMSP['34DMSP_no_enz_deg']/1000)
ax1.plot(df_DMSP['Time'], df_DMSP['32DMSP_no_enz_deg']/1000)
ax1.plot(df_DMSP['Time'], df_DMSP['DMSP_total_no_enz_deg']/1000)
#Set axes labels
ax1.set_xlabel ('Time (min)')
ax1.set_ylabel ('DMSP ($\mu$M)')
#Set axes limits and tick marks
ax1.set_xticks(range(0,70,10))
ax1.set_xlim(0,60)
ax1.set_yticks(range(0,220,30))
ax1.set_ylim(0,210)
#Add legend
ax1.legend(['${34}^{DMSP}$', '${32}^{DMSP}$','DMSP total'])

#Create second plot (right; change in 34DMSP, 32DMSP and total DMSP relative to the fraction of substrate remaining)
ax2.plot(df_DMSP['-ln_fr_no_enz_deg'], df_DMSP['34DMSP_no_enz_deg']/1000)
ax2.plot(df_DMSP['-ln_fr_no_enz_deg'], df_DMSP['32DMSP_no_enz_deg']/1000)
ax2.plot(df_DMSP['-ln_fr_no_enz_deg'], df_DMSP['DMSP_total_no_enz_deg']/1000)
#Set axes labels
ax2.set_xlabel ('-ln $f_{R}$')
ax2.set_ylabel ('DMSP ($\mu$M)')
#Set axes limits and tick marks
ax2.set_xticks(range(0,4))
ax2.set_xlim(0,3.5)
ax2.set_yticks(range(0,220,30))
ax2.set_ylim(0,210)
#Add legend
ax2.legend(['${34}^{DMSP}$', '${32}^{DMSP}$','DMSP total'])
plt.tight_layout(pad=0, w_pad=0.5)

# %%
#Save figure
fig.savefig(f'{homedir}/figures/enz_deg/modelling/change_in_dmsp_without_alma1_deg.pdf', bbox_inches='tight')
# %%
