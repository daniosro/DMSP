# %%
import numpy as np
import pandas as pd
import scipy as sp
import math
import matplotlib.animation as animation
from scipy.integrate import odeint
from numpy import arange
from scipy.integrate import odeint
import scipy.optimize 
from scipy.optimize import leastsq
from math import exp
from collections import OrderedDict
from sklearn.linear_model import LinearRegression
pd.options.mode.chained_assignment = None
import git

# Find home directory for repo
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir

# %% 
# Load data
df_data =  pd.read_csv (f'{homedir}/data/raw/EA/20190917_EA.csv')

# Keep only rows with Peak Nr=2, since Peak Nr=1 
# and Peak Nr=3 correspond to the reference gas
#Select rows to drop
i1 = df_data[df_data['Peak_Nr'] == 1].index
i3 = df_data[df_data['Peak_Nr'] == 3].index
#Drop rows in place
df_data.drop(i1, inplace=True)
df_data.drop(i3, inplace=True)
df_data.head()

# %% 
# Import table with the raw data for blanks
df_blanks = pd.read_csv(f'{homedir}/data/raw/EA/20190917_blanks.csv')

# Remove rows with missing data
df_blanks.dropna(how='any', inplace=True)

# Keep only rows with Peak Nr=2
#Select rows to drop
i1 = df_blanks[df_blanks['Peak_Nr'] == 1].index
i3 = df_blanks[df_blanks['Peak_Nr'] == 3].index
#Drop rows in place
df_blanks.drop(i1, inplace=True)
df_blanks.drop(i3, inplace=True)
df_blanks.head()

# %% 
# Import linearity data
df_lin =  pd.read_csv(f'{homedir}/data/raw/EA/20190917_linearity.csv')
df_lin.head()

# %% 
# Blank correction
# Get the average area all for the blanks
blank_area_all_average = df_blanks["Area_All"].mean()

# Get the average d34S for the blanks
blank_d34s_average = df_blanks["d34S"].mean()

# Append a column to the data dataframe with the 
# correction of area by blank: substraction of the
# area of the measurement by the area of the blank

df_data['Correction_of_area_blank'] = \
    np.subtract(df_data['Area_All'], blank_area_all_average)

# Correction of d34S by blank: 
    
# Get the product of the d34S by the area of each measurement
num1 = np.multiply (df_data['d34S'], df_data['Area_All'])
# Get the product of the d34S by the area of the average of the blanks
num2 = blank_d34s_average * blank_area_all_average  
#Subtract the product of d34S by the product of the blanks
num = np.subtract (num1, num2)
#Divide by the area corrected by blank and
#append column to the data dataframe
df_data['Correction_of_d34S_blank'] = \
    np.divide (num, df_data['Correction_of_area_blank'])
# %% 
# Linearity correction
# Create variables for linear regression
# x will be the amplitude of the peak of mass 64 = ^32S + ^16O + ^16O
# We divided by 1000 to convert volts ot milivolts
x = (df_lin['Ampl_64'].values/1000).reshape((-1,1))
# y will be the values of d34S
y = (df_lin['d34S'])

# Create model of linear regression
model = LinearRegression().fit(x,y)

#Determine R square, intercept and slope
r_sq = model.score(x,y)
intercept = model.intercept_
s = model.coef_
slope=s[0]

# Calculate an amplitude difference by centering around an arbitrary 
# value in the Ampl_64 column of the data dataframe  
num= np.subtract(df_data['Ampl_64'], 700)
ampl_difference = np.divide (num, 1000)

# Calculate the amplitude correction factor
ampl_correction_factor = (slope * ampl_difference)+intercept

# Correct the d34S corrected by amplitude by
# subtracting the blank corrected data by the amplitude correction
# factor and append the column to the data dataframe
df_data['Correction_of_d34S_by_amplitude'] = \
    np.subtract(df_data['Correction_of_d34S_blank'], ampl_correction_factor)

# %% 
#Standard correction

#Create a category column according to the area the peak of each sample and standard
# Correct each sample by standards with similar areas

#Create standard correction group empty list
Std_group = []
# Classify the area of the samples by category 
# Loop through dataframe rows
for row in df_data['Area_All']:
    if row <= 30:
        Std_group.append ('low')
    else:
        Std_group.append ('high')     
# Append to data dataframe 
df_data['Std_group'] = Std_group
df_data.head()

# %%
#Create a dataframe for standards and ditch the outliers
# Create df with only standard data
df_standards = df_data[(df_data.Type == 'Standard')] 
#Sort values by ID and area
df_standards = df_standards.sort_values(['Identifier', 'Area_All'])

#Drop outliers in place
df_standards.drop(df_standards.index[14], inplace=True)
df_standards.drop(df_standards.index[11], inplace=True)
df_standards.drop(df_standards.index[2], inplace=True)

# Calculate the slope and intercept for calculated vs. true value of the stds

# Append true value column to the standard table
# For sulfanilamide
df_standards.loc[df_standards['Identifier'] \
                  == 'Sulfanilamide', 'True_d34S'] = 2.42
# For seawater
df_standards.loc[df_standards['Identifier'] \
                  == 'SW', 'True_d34S'] = 21
                       
# Group by Std_group
df_stdgroup = df_standards.groupby (['Std_group'])

df_standards.head()

# %%
#Linear regression for the standards of each area group
# Get the slope and intercept for the standards of each group
# Define column names
names = ['Std_group', 'R squared', 'Intercept', 'Slope']

# Initialize empty dataframe to save fit results
df_linreg_stds = pd.DataFrame(columns=names)
                 
# Create variables for linear regression
#Loop through standard groups                       
for group, data in enumerate (df_stdgroup):
    #x will be the values of d34S of the stds corrected by blanks and linearity
    x_std = data[1].Correction_of_d34S_by_amplitude.values.reshape((-1,1))
    #y will be the true values of d34S of each standard
    y_std = data[1].True_d34S
    # Create model
    model = LinearRegression().fit(x_std,y_std)
    #Determine R square, intercept and slope
    r_sq_stds = model.score(x_std,y_std)
    intercept_stds = model.intercept_
    s1 = model.coef_
    slope_stds = s1[0]
    # Store parameters and group as list
    params = (data[1].Std_group.unique(), r_sq_stds, intercept_stds, slope_stds)
 
    # Convert list to pandas Series
    series = pd.Series(params, index=names)   
    # Append parameters to dataframe
    df_linreg_stds = df_linreg_stds.append(series, ignore_index=True)
#Round the values of the dataframe to two decimal places
df_linreg_stds = df_linreg_stds.round(2)     
df_linreg_stds 

# %%
#Apply the standard correction to the samples based on their area
#Apply corrections by true value and area

# Initialize lists to save values
slopes = []
intercepts = []
Correction_of_d34S_by_true_value = []

#loop through rows in dataframe
for index, row in df_data.iterrows():
    # Extract standard group
    Std_group = row.Std_group
    # Extract slope and intercept
    slope = df_linreg_stds[df_linreg_stds.Std_group == Std_group].Slope.iloc[0]
    intercept = df_linreg_stds[df_linreg_stds.Std_group == Std_group].Intercept.iloc[0]
    slopes.append(slope)
    intercepts.append(intercept)
    # Compute corrected concentration
    Correction_of_d34S_by_true_value.append(intercept + slope * row.Correction_of_d34S_by_amplitude)
    
# Append values to dataframe
df_data['Correction_of_d34S_by_true_value'] = Correction_of_d34S_by_true_value

df_data.head()

# %%
# Create sample dataframe
df_samples = df_data[(df_data.Type == 'Sample')]
df_samples.head()

# %%
# Export data table
df_samples.to_csv(f'{homedir}/data/processed/EA/20190917_EA.csv')

# %%
# Calculate the analytical repeatibility of the measurements
# Update standard df
# Create a dataframe for only standards from the updated data dataframe
df_standards = df_data[(df_data.Type == 'Standard')]
#Sort values
df_standards = df_standards.sort_values(['Identifier', 'Area_All'])

#Group standards by identifier and amount
grouped_standards = df_standards.groupby(['Identifier', 'Amount'])

# Determine the mean of each standard  and rename the series
mean_stds = grouped_standards['Correction_of_d34S_by_true_value'].mean()
mean_stds = mean_stds.rename("d34S_mean")

# Determine the standard deviation of each standard  and rename the series
std_dev_stds = grouped_standards['Correction_of_d34S_by_true_value'].std()
std_dev_stds = std_dev_stds.rename("d34S_stdev")

#Pass series to individual dataframes
df_mean_stds=mean_stds.to_frame()
df_std_dev_stds=std_dev_stds.to_frame()

#merge the mean and standard deviation dataframes
df_anrep = pd.merge(df_mean_stds, df_std_dev_stds,  how='outer', on=['Identifier', 'Amount'])

#Reset index of the dataframe
df_anrep = df_anrep.reset_index()

#Add column of true value of the standards
df_anrep.loc[df_anrep['Identifier'] \
                  == 'Sulfanilamide', 'True_d34S'] = 2.42
df_anrep.loc[df_anrep['Identifier'] \
                  == 'SW', 'True_d34S'] = 21

#Determine the accuracy by subtracting the true value from the average value of each standard
df_anrep ['Accuracy'] = abs(df_anrep ['True_d34S'] - df_anrep ['d34S_mean'])
df_anrep

# %%
