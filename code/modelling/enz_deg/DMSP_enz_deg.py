# %%
#For numerical calculations
import numpy as np
from scipy.integrate import odeint
import operator 
from scipy.optimize import fsolve
from scipy.stats import linregress
import matplotlib.pyplot as plt
import pandas as pd
import git

# Find home directory for repo
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir

# %% 
#Let's define the function isotopes, which calculates the amount of 
# 34DMSP and 32DMSP for a given DMSP concentration.

def isotopes(c, r_34_std = 0.0450045, delta_in = 14.3):
    
    '''
    Function that splits a DMSP concentration in the proportion of
    heavy and light atoms.
    Parameters
    ----------
    c: float.
    Concentration of DMSP.
    r_34_std: float.
    34R of the VCDT standard. Error = +-0.0093
    delta_in: float.
    Delta 34S of the DMSP used in the assays for the assays, from Sigma.
    Returns
    -------
    Concentration of 34DMSP, concentration of 32DMSP
    '''
        
    # Calculate 34R in DMSP from Eq. 3
    r_34_dmsp =r_34_std *(1+(delta_in/1000))
    # Calculate 32_DMSP from Eq. 2
    dmsp_32 = c/(1 + r_34_dmsp)
    # Calculate 34_DMSP from Eq. 1
    dmsp_34 = c - dmsp_32
        
    return [dmsp_34,dmsp_32]

# %% 
# Now, let's define a function that will numerically integrate the change in the 
# concentrations of 34DMSP and 32DMSP over time.

def dmsp_enz_deg(
    c, 
    t,   
    alpha,
    vmax,
    vmax_32,
    kappa_32,
    k
):
    """
    Function that computes dD32_dt and dD34_dt of DMSP
    Parameters
    ----------
    c: float.
    Concentration of DMSP in nM.
    t: int
    Integration time in min.
    alpha: float.
    Alpha for cleavage by DddP from this study. 
    vmax: float.
    Vmax for cleavage by DddP, calculated from the K M that the enzyme should have to 
    exhibit the pattern of d34S DMSP vs. time, in nM/min/nM enzyme
    Vmax_d: float.
    km: float.
    K M that the enzyme should have to exhibit the pattern of d34S DMSP vs. time, in nM.
    k: float.
    Degradation rate of the enzyme, in min^-1.

    Returns
    -------
    The dD32_dt and dD34_dt of DMSP
    """

    # Unpack isotopes
    enzyme, dmsp_34, dmsp_32 = c
    
    #Calculate vmax_34 assuming that Vmax total = Vmax_32 + Vmax_34
    #This assumption would only hold true at saturation
    vmax_34 = vmax-vmax_32
    
    #Determination of kappa 32 from kappa 34 and the fractionation factor
    kappa_34 = kappa_32 * alpha
    
    # Calculate dD34_dt
    dD34_dt = - ((kappa_34 * enzyme * (vmax_34 * enzyme * dmsp_34/((vmax_34 * enzyme)+(kappa_34 * enzyme * dmsp_34)))))
    
    # Calculate dD32_dt
    dD32_dt = - ((kappa_32 * enzyme * (vmax_32 * enzyme * dmsp_32/((vmax_32 * enzyme)+(kappa_32 * enzyme * dmsp_32)))))
    
    #Calculate dE_dt
    dE_dt = -k*enzyme

    return [dE_dt, dD34_dt, dD32_dt]

# %% 
#Let's define the parameters for integration
# Define parameters for integration

# For DddP (this study)
alpha = (-3.97/1000)+1
# Initial concentration of DMSP in nM
dmsp_init = 207*1000
#Get the isotopic composition of initial DMSP
dmsp_iso_init = isotopes(dmsp_init)
#c[1] is the starting enzyme concentration in nM and c[2] is the starting DMSP concentration in nM
#For DddP
c = (3.479,*dmsp_iso_init)
#For DddP (in nM)
#This is the value of K_M that allows to better reproduce the 
# slope with coherent d34S values of remaining DMSP. Units are nM.
km = 2E9
#Vmax in nM/min/nM enzyme. It will be calculated by reorganizing Eq. 2
# d[DMSP]/dt was assumed to be equal to ([second point] - [first point])/min/enzyme concentration
# In this case, the second point was [DMSP] at t=0 and the first point was [DMSP] at t = 3 min.
vmax=(17000*(dmsp_init+km))/(dmsp_init*c[0])
# Define the 34R of the standard (VCDT)
r_34_std = 0.0450045
#Define the time points for integration
t = np.linspace(0, 53, 50)
#This is the value of V_max calculated using the K_M above, in nM/min/nM enzyme
vmax
# %% 
#This is the ratio of Vmax to K_M (kappa)
kappa = vmax/km
kappa
# %% 
#Now, let's integrate the equations above
#Determine the fractionation by the enzyme when k=0

#k is the rate constant of DMSP loss in min^-1
k = 0
# This is the minimum value that Vmax_32 can have to explain our results,
# considering that 32S is much more abundant than 34S, so Vmax ~ 32Vmax
vmax_32=vmax*0.8
#We will assume that kappa_32 = kappa
kappa_32=kappa
# Group all parameters into args
args = (
    alpha,
    vmax,
    vmax_32,
    kappa_32,
    k
)
# Integrate numerically with ODEint
dmsp_iso = odeint(dmsp_enz_deg, c, t, args=args)

# Split the values returned by the integration, using k=0
enzyme_change_0= dmsp_iso[:,0]
DMSP_34_0 = dmsp_iso[:,1]
DMSP_32_0 = dmsp_iso[:,2]

# Calculate 34R for k=0
r_34_0 = DMSP_34_0/DMSP_32_0

# Calculate d34S and total DMSP for k=0
d34s_0 = (np.divide(r_34_0,r_34_std)-1) * 1000
total_0 = DMSP_34_0+DMSP_32_0

# Calculate the fraction of substrate remaining for k=0
f_r_0 = total_0/dmsp_init
#Determine -ln (f_R) for k=0
minusnatlog_f_r_0 = -np.log(f_r_0)
#Get the slope of the line in a plot of -ln (f_R) vs. the d34s of DMSP,
#equivalent to the fractionation factor without loss of enzyme activity
slope_0, _, _, _, _ = linregress(minusnatlog_f_r_0, d34s_0)

# %% 
#Determine the fractionation by the enzyme when k is different from 0
#k is the rate constant of DMSP loss in min^-1
k =0.08

# Group all parameters into args
args = (
    alpha,
    vmax,
    vmax_32,
    kappa_32,
    k
)

# Integrate numerically with ODEint
dmsp_iso = odeint(dmsp_enz_deg, c, t, args=args)

# Split the values returned by the integration, using k=!0
enzyme_change_k= dmsp_iso[:,0]
DMSP_34_k = dmsp_iso[:,1]
DMSP_32_k = dmsp_iso[:,2]

# Calculate 34R for k=!0
r_34_k = DMSP_34_k/DMSP_32_k

# Calculate d34S and total DMSP for k=!0
d34s_k = (np.divide(r_34_k,r_34_std)-1) * 1000
total_k = DMSP_34_k+DMSP_32_k

# Calculate the fraction of substrate remaining for k=!0
f_r_k = total_k/dmsp_init
#Determine -ln (f_R) for k=!0
minusnatlog_f_r_k = -np.log(f_r_k)
#Get the slope of the line in a plot of -ln (f_R) vs. the d34s of DMSP,
#equivalent to the fractionation factor with loss of enzyme activity
slope_k, _, _, _, _ = linregress(minusnatlog_f_r_k, d34s_k)

# %% 
#Let's create a dataframe to store the variables from the inegration before that we
# are interested in plotting
#Let's create a dataframe to store the variables from the inegration before that we
# are interested in plotting
df_DMSP = pd.DataFrame(list(zip(t, DMSP_34_k, DMSP_32_k, total_k,
                               minusnatlog_f_r_k, d34s_k, DMSP_34_0, DMSP_32_0, 
                               total_0, minusnatlog_f_r_0, d34s_0)), 
               columns =['Time', '34DMSP_enz_deg', '32DMSP_enz_deg',
                        'DMSP_total_enz_deg','-ln_fr_enz_deg', 'd34s_enz_deg',
                         '34DMSP_no_enz_deg','32DMSP_no_enz_deg',
                         'DMSP_total_no_enz_deg','-ln_fr_no_enz_deg','d34s_no_enz_deg']) 
df_DMSP.head()
# %% 
#Let's export the previous dataframe
df_DMSP.to_csv(f'{homedir}/data/modelling/DMSP_enz_deg.csv')

# %%
