# %%
import os
import numpy as np
import scipy
import scipy.optimize
import pandas as pd
import git
from collections import OrderedDict

# Find home directory for repo
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir

# Import plotting features
import matplotlib.pyplot as plt
import seaborn as sns
import bokeh.plotting
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, FactorRange
from bokeh.palettes import all_palettes
from bokeh.io import export_png
# %%
# Set plot style
sns.set_style("ticks")
sns.set_palette("colorblind", color_codes=True)
sns.set_context("paper")

# load data
df_enzymes = pd.read_csv(f'{homedir}/data/processed/genetics/stab8_cursonetal_2018_tidy.csv')
df_enzymes.head()

# %%
# Make stacked bar plot for marine pelagic biomes following 
# the Longhurst code: https://en.wikipedia.org/wiki/Longhurst_code.

# Create a dataframe with the average of transcripts for each marine pelagic biome.
# Get averages for the total number of transcripts per million sequences for each marine pelagic biome and enzyme
df_averages2_ttpms = df_enzymes.groupby(['Marine_pelagic_biomes _Longhurst_2007', 'Enzyme']).Total_transcripts_per_million_sequences.agg(['mean','std']).reset_index()
# Pivot table to classify the means by biome and enzyme
df_averages2_ttpms = df_averages2_ttpms.pivot_table(index=['Marine_pelagic_biomes _Longhurst_2007'], columns='Enzyme', values='mean').reset_index()
df_averages2_ttpms
# %%
#Make stacked plot for enzymes with marine pelagic biome
#Define data source
source = ColumnDataSource(data=df_averages2_ttpms)

#Define lists of unique values
marine_pelagic_biomes = df_enzymes[('Marine_pelagic_biomes _Longhurst_2007')].unique()
enzymes = df_enzymes[('Enzyme')].unique()
enzyme_names = enzymes.tolist()

#Define color palette
color = all_palettes['Colorblind'][len(enzymes)]

# Make figure
p = figure(x_range=marine_pelagic_biomes,
           width=800,
           height=500,
           x_axis_label='Marine pelagic biomes',
           y_axis_label='Total transcripts per million sequences',     
           tools="")

# Add stacked bars
p.vbar_stack(enzymes, x='Marine_pelagic_biomes _Longhurst_2007', 
            width=0.9, source = source, fill_color=color,
             legend_label=enzyme_names)


#Styling
p.xaxis.axis_label_text_font_size = "16pt"
p.yaxis.axis_label_text_font_size = "16pt"
p.xaxis.major_label_text_font_size = "12pt"
p.yaxis.major_label_text_font_size = "14pt"

bokeh.io.show(p)

#Export figure
export_png(p, filename=f'{homedir}/figures/genetics/metatranscriptomics/stab8_curson2018_bact_biomes.png')

# %%
