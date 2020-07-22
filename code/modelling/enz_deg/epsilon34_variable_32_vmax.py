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
# Now, let's define a function that will numerically integrate the change 
# in the concentrations of 34DMSP and 32DMSP over time.

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
    Vmax for cleavage by DddP, calculated from the K M that the enzyme 
    should have to exhibit the pattern of d34S DMSP vs. time, in nM/min/nM enzyme
    Vmax_d: float.
    km: float.
    K M that the enzyme should have to exhibit the pattern of d34S DMSP 
    vs. time, in nM.
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
    dD34_dt = - ((kappa_34 * enzyme * (vmax_34 * enzyme * dmsp_34/
    ((vmax_34 * enzyme)+(kappa_34 * enzyme * dmsp_34)))))
    
    # Calculate dD32_dt
    dD32_dt = - ((kappa_32 * enzyme * (vmax_32 * enzyme * dmsp_32/
    ((vmax_32 * enzyme)+(kappa_32 * enzyme * dmsp_32)))))
    
    #Calculate dE_dt
    dE_dt = -k*enzyme

    return [dE_dt, dD34_dt, dD32_dt]
# %%
# Define parameters for integration, such that the only thing that varies is
# 32_Vmax

# For DddP (this study)
alpha = (-3.97/1000)+1
# Initial concentration of DMSP in nM
dmsp_init = 207*1000
#Get the isotopic composition of initial DMSP
dmsp_iso_init = isotopes(dmsp_init)
#c[1] is the starting enzyme concentration in nM and c[2] is the starting DMSP 
# concentration in nM
#For DddP
c = (3.479,*dmsp_iso_init)
#For DddP (in nM)
#This is the value of K_M that allows to better reproduce the 
# slope with coherent d34S values of remaining DMSP. Units are nM.
km = 2E9
#Vmax in nM/min/nM enzyme. It will be calculated by reorganizing Eq. 2
# d[DMSP]/dt was assumed to be equal to ([second point] - 
# [first point])/min/enzyme concentration
# In this case, the second point was [DMSP] at t=0 and the first point 
# was [DMSP] at t = 3 min.
vmax=(17000*(dmsp_init+km))/(dmsp_init*c[0])
# Define the 34R of the standard (VCDT)
r_34_std = 0.0450045
#Define the time points for integration
t = np.linspace(0, 53, 50)
#Create an array with a range of enzyme degradation rates
k = np.linspace (0,1,100)

#Define Vmax and 32^Vmax for cleavage by DddP in nmol/min/mg enz
vm=vmax
vmax_32=vmax*0.8
#Define K_M
km = 2E9
# Define the ratio of Vmax to K_M (kappa)
kappa = vmax/km
#Set 32_kappa equal to kappa
kappa_32=kappa

#Create array of 32^Vmax for cleavage by DddP in nmol/min/mg enz
vmax_32 = np.linspace (1E6,1E9,100)
#Define Vmax in nmol/min/mg enz
vm=2E9
#Define enzyme degradation rate in min^-1
k=0.08

# %%
#Let's perform the integration
#Initialize empty lists
slopes = []
vmax_32_x = []
#Loop through the enzyme degradation rates
for i in vmax_32:
    # Group all parameters into args
    args = (
    alpha,
    vm,
    i,
    kappa_32,
    k
    )
    # Integrate numerically with ODEint
    dmsp_iso = odeint(dmsp_enz_deg, c, t, args=args)
    # Split the values returned by the integration
    enzyme_change= dmsp_iso[:,0]
    DMSP_34 = dmsp_iso[:,1]
    DMSP_32 = dmsp_iso[:,2]
    # Calculate 34R
    r_34 = DMSP_34/DMSP_32
    # Calculate d34S and total DMSP
    d34s = (np.divide(r_34,r_34_std)-1) * 1000
    lin_app_d34s = 1000*np.log(1+d34s/1000)
    total = DMSP_34+DMSP_32
    #Calculate the fraction of reactant remaining
    f_r = total/dmsp_init
    #Determine -ln (f_R)
    minusnatlog_f_r = -np.log(f_r)
    
    #Get the slope of the line in a plot of -ln (f_R) vs. the linear 
    # approximation to the d34s of DMSP,
    #equivalent to the fractionation factor
    slope, _, _, _, _ = linregress(minusnatlog_f_r, lin_app_d34s)
    #Append the fractionation factors for different Vmax to the previously 
    # created empty list
    slopes.append(slope)
    #Append the Vmax to the previously created empty list
    vmax_32_x.append (i)
    
#Send slopes and 32^V_max to a new dataframe
df_34e_variable_32vmax = pd.DataFrame(list(zip(vmax_32_x,slopes)), 
               columns =['32vmax','34epsilon']) 
df_34e_variable_32vmax.head()
# %%
#Export dataframe to csv
df_34e_variable_32vmax.to_csv(
    f'{homedir}/data/modelling/34e_variable_32vmax.csv')