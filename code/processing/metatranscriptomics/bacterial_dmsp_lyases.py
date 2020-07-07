# %%
import os
import numpy as np
import scipy
import scipy.optimize
import pandas as pd
from collections import OrderedDict

# Import plotting features
import matplotlib.pyplot as plt
import seaborn as sns

# %%
# Set plot style
sns.set_style("ticks")
sns.set_palette("colorblind", color_codes=True)
sns.set_context("paper")

# load data
df_stab8_cursonetal_2018 = pd.read_csv('/Users/daniosro/git/DMSP/data'
'/stab8_cursonetal_2018.csv')
df_stab8_cursonetal_2018.head()

# %%
