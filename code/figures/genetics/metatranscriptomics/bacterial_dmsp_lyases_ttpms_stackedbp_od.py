# %%
import numpy as np
import pandas as pd
import git

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
# Make stacked bar plot for ocean depths.
# DCM = Depth of chlorophyll maximum, 
# MES = mesopelagic ocean, MIX = mixed layer depth, 
# SRF = surface ocean. 

# Create a dataframe with the average of transcripts for each depth region.
# Get averages for the total number of transcripts per million sequences for each enzyme and depth zone
df_averages_ttpms = df_enzymes.groupby(['Enzyme', 'Depth_zone']).Total_transcripts_per_million_sequences.agg(['mean','std']).reset_index()
# Pivot table to classify the means by depth and enzyme
df_averages_ttpms = df_averages_ttpms.pivot_table(index=['Enzyme'], columns='Depth_zone', values='mean').reset_index()

df_averages_ttpms
# %%
#Make stacked plot for enzymes with depth zones
#Define data source
source = ColumnDataSource(data=df_averages_ttpms)

#Define lists of unique values
names = df_enzymes['Depth_zone'].unique().tolist()
depth_zones = df_enzymes[('Depth_zone')].unique()
enzymes = df_enzymes[('Enzyme')].unique()

#Define color palette
color = all_palettes['Colorblind'][len(depth_zones)]

# Make figure
p = figure(x_range=enzymes,
           width=800,
           height=500,
           x_axis_label='Gene',
           y_axis_label='Total transcripts per million sequences',     
           tools="")

# Add stacked bars
p.vbar_stack(depth_zones, x='Enzyme', 
            width=0.9, source = source, fill_color=color,
             legend_label=names)

#Styling
p.xaxis.axis_label_text_font_size = "16pt"
p.yaxis.axis_label_text_font_size = "16pt"
p.xaxis.major_label_text_font_size = "14pt"
p.yaxis.major_label_text_font_size = "14pt"

# Eliminate grid
p.xgrid.grid_line_color = None
p.ygrid.grid_line_color = None
#Set outline
p.outline_line_width = 1
p.outline_line_color = "black"
#Eliminate minor ticks
p.yaxis.minor_tick_line_color = None
#Set font in axes to normal
p.xaxis.axis_label_text_font_style = "normal"
p.yaxis.axis_label_text_font_style = "normal"


bokeh.io.show(p)
#Export figure
export_png(p, filename=f'{homedir}/figures/genetics/metatranscriptomics/stab8_curson2018_bact_depths.png')

# %%
