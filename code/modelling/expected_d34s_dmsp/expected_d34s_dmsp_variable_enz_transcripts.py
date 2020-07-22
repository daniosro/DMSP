# %% 
#For numerical calculations
import numpy as np
import operator
from scipy.integrate import odeint 
from scipy.optimize import fsolve
import pandas as pd
pd.set_option('display.precision',5)
import git

# Find home directory for repo
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir
# %% 
#Let's start by defining a function that determines the isotopic composition
#of DMSP for a particular concentration or flux of DMSP
def isotopes(c, r_34_std = 0.0450045, delta_in = 17):
    '''
    Function that splits a DMSP concentration in the proportion of
    heavy (34) and light (32) atoms.
    Parameters
    ----------
    c: float.
    Concentration of DMSP in nM.
    r_34_std: float.
    34R of the VCDT standard. Error = +-9.3
    delta_in: float.
    Delta 34S of newly synthesized DMSP. Assumed.
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
#Now, we will define a function that calculates the change in 32DMSP and 34DMSP
#over time, when DMSP is degraded by DmdA (DMSP demethylase), Alma1 (eukaryotic 
# DMSP lyase) and DddP (most abundant bacterial DMSP lyase)
def dmsp_system_3_comp_mm(
    c, 
    t,
    f_total_in, 
    alpha_d,
    alpha_c1,
    alpha_cp,
    kappa_32_d, 
    kappa_32_c1,
    kappa_32_cp,
    vmax_34_d,
    vmax_34_c1,
    vmax_34_cp,
    vmax_32_d,
    vmax_32_c1,
    vmax_32_cp,
    transcripts_d,
    transcripts_c1,
    transcripts_cp
):
    """
    Function that computes dD32_dt and dD34_dt of DMSP
    Parameters
    ----------
    c: float.
    Concentration of DMSP in nM.
    t: int
    Integration time in min.
    f_total_in: float. 
    [DMSP] that enters the system in nmol/min. Range: 0.2-2.5 nmol/l/h. From Simo et al. (2009)
    alpha_d: float.
    Alpha for demethylation from this study. 
    Epsilon has an error of +-0.17
    alpha_c1: float.
    Alpha for cleavage by Alma1 from this study. 
    Epsilon has an error of +-0.06
    alpha_cp: float.
    Alpha for cleavage from this study. 
    Epsilon has an error of +-0.14 
    kappa_32_d:float.
    Ratio of 32^V_max to 32^K_M for demethylation in l/min/mg enz. 
    kappa_32_c1:float.
    Ratio of 32^V_max to 32^K_M for cleavage by Alma1 in l/min/mg enz. 
    kappa_32_cp:float.
    Ratio of 32^V_max to 32^K_M for cleavage by DddP in l/min/mg enz. 
    vmax_34_d: int.
    34^V_max for demethylation in nmol/min/mg enz.
    vmax_34_c1: int.
    34^V_max for cleavage by Alma1 in nmol/min/mg enz.
    vmax_34_cp: int.
    34^V_max for cleavage by DddP in nmol/min/mg enz.
    vmax_32_d: int.
    32^V_max for demethylation in nmol/min/mg enz.
    vmax_32_c1: int.
    32^V_max for cleavage by Alma1 in nmol/min/mg enz.
    vmax_32_cp: int.
    32^V_max for cleavage by DddP in nmol/min/mg enz.
    transcripts_d: int.
    Concentration of DmdA (demethylase) in the ocean from Varaljay et al. (2015).
    Units of mRNA/l.
    transcripts_c1: int.
    Concentration of DMSP lyase in the ocean from approximation based on Vorobev et al. (2020)
    Units of mRNA/l.
    transcripts_cp: int.
    Concentration of DMSP lyase in the ocean from Varaljay et al. (2015)
    Units of mRNA/l.

    Returns
    -------
    The dD32_dt and dD34_dt of DMSP
    """
    # Unpack isotopes
    dmsp_34, dmsp_32 = c
    #Flux in of each DMSP
    f_34_in, f_32_in = isotopes(f_total_in)
    
    #Determination of kappa 34 for cleavage by Alma1 from kappa 32 and the fractionation factor
    kappa_34_c1 = kappa_32_c1 * alpha_c1 
    #Determination of kappa 34 for cleavage by DddP from kappa 32 and the fractionation factor
    kappa_34_cp = kappa_32_cp * alpha_cp 
    #Determination of kappa 34 for demethylation from kappa 32 and the fractionation factor
    kappa_34_d = kappa_32_d * alpha_d 
    
    # Calculate enzyme concentrations in units of mg/L
    # Final 1000 is the conversion of grams to milligrams. 

    dmda = prots_mrna*transcripts_d*(1/avog_n)*dmda_w*1000
    alma1 = prots_mrna*transcripts_c1*(1/avog_n)*alma1_w*1000
    dddp = prots_mrna*transcripts_cp*(1/avog_n)*dddp_w*1000
    
    #Calculate flux out of 32^DMSP in units of nmol/min 
    f_32_out = ((kappa_32_c1 * alma1 * (vmax_32_c1 * dmsp_32/((vmax_32_c1)+(kappa_32_c1 * dmsp_32))))+
                (kappa_32_cp * dddp * (vmax_32_cp * dmsp_32/((vmax_32_cp)+(kappa_32_cp * dmsp_32))))+
                (kappa_32_d * dmda * (vmax_32_d * dmsp_32/((vmax_32_d)+(kappa_32_d* dmsp_32)))))*ocean_vol
    
    # Calculate dD32_dt in units of nmol/l/min
    dD32_dt = (f_32_in - f_32_out)/ocean_vol
    
    #Calculate flux out of 34^DMSP in units of nmol/min
    f_34_out = ((kappa_34_c1 * alma1 * (vmax_34_c1 * dmsp_34/((vmax_34_c1)+(kappa_34_c1 * dmsp_34))))+
                (kappa_34_cp * dddp * (vmax_34_cp * dmsp_34/((vmax_34_cp)+(kappa_34_cp * dmsp_34))))+
                (kappa_34_d * dmda * (vmax_34_d * dmsp_34/((vmax_34_d)+(kappa_34_d * dmsp_34)))))*ocean_vol
    # Calculate dD34_dt in units of nmol/l/min
    dD34_dt = (f_34_in - f_34_out)/ocean_vol
    
    return [dD34_dt, dD32_dt]

# %%
#Let's define important variables that will be incorporated in the calculation, 
# as well as integration parameters.
#Number of proteins/mRNA. 
prots_mrna = 1000
#Avogadro number
avog_n=6.022E23

#Weight of DmdA in g/mol
dmda_w=39.5E3
#Weight of DddP in g/mol
dddp_w=50E3
#Weight of Alma1 in g/mol
alma1_w=160E3

#Other variables
#Ocean volume in liters
ocean_vol = 1.3E21
#Order of magnitude estimate of the proportion of enzyme relative to total protein 
# in the cell
total_prot_enzyme_ratio = 1/1E-6
#Average DMSP concentration in nM
dmsp_ocean = 10
#Isotopic composition of DMSP in the ocean from the isotopes function
isotopes_ini=isotopes(dmsp_ocean)

# Set parameters for integration
#Initial concentration of DMSP in nM
c = isotopes(10)
# Flux of DMSP into the ocean in nmol/l/min from Simó et al. (2019)
#(Annual average for the Mediterranean)
f_in = (0.000015/60)
# Volume of the surface ocean in liters, assuming a depth of 200 m
surface_ocean_vol = 7.24E19
# Flux of DMSP into the ocean in nmol/l/min
f_total_in = f_in*surface_ocean_vol

#Fractionation factors determined in this study
#For DmdA
alpha_d = (-2.72/1000)+1
# For Alma1
alpha_c1 = (-1.18/1000)+1
# For DddP
alpha_cp = (-3.97/1000)+1

# Define the 34R of the standard
r_34_std = 0.0450045
#Define the time points for integration
t = np.linspace(0, 1E10, 1000)

# %%
# Based on the lower estimates of the kinetic parameters for the enzymes in the 
# native organisms reported by Jonkers et al. (2000) and Stefels et al. (2007), 
# we will calculate the lower limit for V_max and V_max/ K_M of each 
# enzyme for open ocean conditions
#V_max in units of nmol/min/mg cell protein
# K_M in nM
#DmdA
vmax_d = 0.089*1000 
km_d = 15500
#Alma1
vmax_c1 = 0.011*1000 
km_c1 = 500000
#DddP
vmax_cp = 0.0144*1000 
km_cp = 70000

#Calculation of Vmax in nmol/min/mg enzyme
vmax_enz_c1 = (vmax_c1*total_prot_enzyme_ratio)
vmax_enz_d = (vmax_d*total_prot_enzyme_ratio)
vmax_enz_cp = (vmax_cp*total_prot_enzyme_ratio)

#Determination of the real ratio of Vmax/kM in l/min/mg enz
ratio_cp = vmax_enz_cp/(km_cp)
ratio_c1 = vmax_enz_c1/(km_c1)
ratio_d = vmax_enz_d/(km_d)

#Append to dataframe
d_low = {'Enzyme': ['Alma1', 'DddP', 'DmdA'],
     '$V_{max} (nmol/min/mg enz)$': [vmax_enz_c1, vmax_enz_cp, vmax_enz_d], 
    '$V_{max}/K_M (l/min/mg\ enz)$': [ratio_c1, ratio_cp, ratio_d]}

df_params_low = pd.DataFrame(data=d_low)
df_params_low

# %%
# Now, we will use the upper estimates of the kinetic parameters for the enzymes in 
# the native organisms upper limit for V_max and V_max/ K_M of each enzyme 
# for open ocean conditions.
#V_max in units of nmol/min/mg cell protein
# K_M in nM
#DmdA
vmax_d = 0.3*1000
km_d = 4100 
#Alma1
vmax_c1 = 0.0835*1000
km_c1 = 12000
#DddP
vmax_cp = 0.0201*1000
km_cp = 20000 

#Calculation of Vmax in nmol/min/mg enzyme
vmax_enz_c1 = (vmax_c1*total_prot_enzyme_ratio)
vmax_enz_d = (vmax_d*total_prot_enzyme_ratio)
vmax_enz_cp = (vmax_cp*total_prot_enzyme_ratio)

#Determination of the real ratio of Vmax/kM in l/min/mg enz
ratio_cp = vmax_enz_cp/(km_cp)
ratio_c1 = vmax_enz_c1/(km_c1)
ratio_d = vmax_enz_d/(km_d)

#Append to dataframe
d_high = {'Enzyme': ['Alma1', 'DddP', 'DmdA'],
     '$V_{max} (nmol/min/mg enz)$': [vmax_enz_c1, vmax_enz_cp, vmax_enz_d], 
    '$V_{max}/K_M (l/min/mg\ enz)$': [ratio_c1, ratio_cp, ratio_d]}

df_params_high = pd.DataFrame(data=d_high)
df_params_high
# %%
# We will perform a numerical integration with ODEINT until the diff. equations that track the 
# change of 32DMSP and 34DMSP over time reach a steady-state.
# We will do so by assuming that the DMSP degrading enzymes have a variable number of transcripts.
# Let's begin by defining the Michaelis-Menten parameters.

#Kappa in l/min/mg enz. Intermediate value
# We will assume that kappa_32 ~ kappa
#DmdA
kappa_32_d=10000
#Alma1
kappa_32_c1 = 2000
#DddP
kappa_32_cp= 600

#V max in nmol/min/mg enz. The values were selected at an intermediate value 
#between the reported extreme values of Vmax for each enzyme. 
# We will assume that in all cases 34^V_max ~ 32^V_max
#Similar results were obtained when V_max ~ 32^V_max
#DmdA
vmax_34_d = 1E7/2
#Alma1
vmax_34_c1 = 4E7/2
#DddP
vmax_34_cp = 9E7/2
#DmdA
vmax_32_d = 1E7/2
#Alma1
vmax_32_c1 = 4E7/2
#DddP
vmax_32_cp = 9E7/2

#Transcripts in mRNA/L.
#Alma1
transcripts_c1 = np.linspace(0,3E7,20)
#DddP
transcripts_cp = np.linspace(0,3E7,20)

# %%
# Perform numerical integration with ODEint
#Create empty lists
DMSP_34 = []
DMSP_32 = []
#Fraction of DMSP degraded by Alma1
fr_alma1=[]
#Fraction of DMSP degraded by DddP
fr_dddp=[]
#Fraction of DMSP degraded by DmdA
fr_dmda = []
#Ttotal DMSP degraded
total = []
#Loop through values of fraction of DMSP cleaved by Alma1
for i in transcripts_c1:
    #Loop through values of fraction of DMSP cleaved by DddP
    for j in transcripts_cp:
        #Calculate the fraction of degradation of DmdA
        transcripts_d=(3E7)-i-j
        #Determine the sum of fraction of DMSP degraded by Alma1 and DddP
        total_temp = i+j+transcripts_d
        #If the previously calculated sum is less than 3E7
        if (total_temp <=3E7 and transcripts_d>=0): 
        #Integrate
            dmsp_iso = odeint(dmsp_system_3_comp_mm, c, t, args= (
                f_total_in,
                alpha_d,
                alpha_c1,
                alpha_cp,
                kappa_32_d, 
                kappa_32_c1,
                kappa_32_cp,
                vmax_34_d,
                vmax_34_c1,
                vmax_34_cp,
                vmax_32_d,
                vmax_32_c1,
                vmax_32_cp,
                transcripts_d,
                i,
                j, 
                )
             )
            #Calculate the fraction of DMSP degraded by each enzyme
            #Alma1
            fr_1 = i/3E7
            #DddP
            fr_p = j/3E7
            #DmdA
            fr_a = transcripts_d/3E7
            #Append 34^DMSP and 32^DMSP to empty lists
            DMSP_34.append(dmsp_iso[-1,0])
            DMSP_32.append(dmsp_iso[-1,1])
            #Append total and fractions degraded by each enzyme in each iteration to empty list
            total.append(total_temp)
            fr_alma1.append(fr_1)
            fr_dddp.append(fr_p)
            fr_dmda.append(fr_a)
# %%
#Let's determine the expected delta 34^s values of DMSP when it is degraded by different fractions of
#each enzyme
# Calculate 34R
r_34 = [x/y for x, y in zip (DMSP_34,DMSP_32)]

#Calculate d34S
#Create empty list
d34s_int = []
#Loop through 34^R
for i in r_34:
    #Calculate d34S
    d34s = ((i/r_34_std)-1)*1000
    #Append lin. app. to empty list
    d34s_int.append(d34s)

#Create dataframe with the information from the integration
df_dmsp_system_transcripts = pd.DataFrame(
    {'fr_alma1': fr_alma1,
     'fr_dddp': fr_dddp,
     'fr_dmda': fr_dmda,
     'd34S': d34s_int
    })
df_dmsp_system_transcripts.head()
# %%
#Export dataframe
df_dmsp_system_transcripts.to_csv(f'{homedir}/data/modelling/expected_d34s_DMSP_variable_enz_transcripts.csv')
