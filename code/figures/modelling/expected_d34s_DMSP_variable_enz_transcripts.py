# %% 
import git
import pandas as pd
import numpy as np

# Find home directory for repo
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir

#Import plotting features
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.tri as tri
import seaborn as sns

# Set plot style
sns.set_style("ticks")
sns.set_palette("colorblind", color_codes=True)
sns.set_context("paper")

# %%
# Load data
df_dmsp_system_transcripts = pd.read_csv(f'{homedir}/data/modelling/expected_d34s_DMSP_variable_enz_transcripts.csv')
df_dmsp_system_transcripts.head()

# %%
#Convert lists to arrays for the plot
fr_dddp = np.asarray(df_dmsp_system_transcripts['fr_dddp'])
fr_dmda = np.asarray(df_dmsp_system_transcripts['fr_dmda'])
fr_alma1 = np.asarray(df_dmsp_system_transcripts['fr_alma1'])
d34s_int = np.asarray(df_dmsp_system_transcripts['d34S'])

# %%
#Let's make the ternary plot for for the fractional contributions of DmdA, DddP and 
# Alma1 to the predicted delta ^34S of DMSP, with a variable number of transcripts of 
# each enzyme.
# Translate the data to Cartesian coordinates
x = 0.5 * ( 2.*fr_dddp+fr_dmda) / (fr_alma1+fr_dddp+fr_dmda )
y = 0.5*np.sqrt(3) * fr_dmda / (fr_alma1+fr_dddp+fr_dmda)

#Create very small intervals for the contours to soften them
tlev = np.linspace(d34s_int.min(),d34s_int.max())

# Create a triangulation out of these points
T = tri.Triangulation(x,y)

# Create figure and define size
fig, ax = plt.subplots(figsize=(2.95, 1.95), dpi=192)

# Plot the contour
plot_tr = plt.tricontourf(x,y,T.triangles,d34s_int, levels=tlev, cmap='viridis')

#Add colorbar
cbar = fig.colorbar(plot_tr,format="%.1f")
#Fix colorbar
cbar_ticks = np.linspace(d34s_int.min(),d34s_int.max(), num=5)
cbar.set_ticks(cbar_ticks)

#Set the right number of ticks on the x axis
ax.set_xticks (np.linspace(0,1,6))

#Set labels 
plt.title('Transcripts/L')
plt.xlabel('$f_{DddP}$')
plt.ylabel('$f_{Alma1\;left} | f_{DmdA\;right}$')
# %%
#Save plot to pdf
fig.savefig(f'{homedir}/figures/modelling/ternary_model_tc.pdf', bbox_inches='tight')