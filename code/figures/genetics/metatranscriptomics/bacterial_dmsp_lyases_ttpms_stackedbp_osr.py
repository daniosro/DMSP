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
# Make stacked bar plot for ocean and sea regions.
# IO = Indian Ocean, NAO = North Atlantic Ocean, 
# NPO = North Pacific Ocean, SAO = South Atlantic Ocean,
# SO = Southern Ocean, SPO = South Pacific Ocean. 

# Create a dataframe with the average of transcripts for each ocean and sea region.
# Get averages for the total number of transcripts per million sequences for each marine pelagic biome and enzyme
df_averages3_ttpms = df_enzymes.groupby(['Ocean_and_sea_regions', 'Enzyme']).Total_transcripts_per_million_sequences.agg(['mean','std']).reset_index()
# Pivot table to classify the means by oceanic region and enzyme
df_averages3_ttpms = df_averages3_ttpms.pivot_table(index=['Ocean_and_sea_regions'], columns='Enzyme', values='mean').reset_index()
df_averages3_ttpms
# %%
#Make stacked plot for enzymes with ocean and sea region
#Define data source
source = ColumnDataSource(data=df_averages3_ttpms)

#Define lists of unique values
ocean_and_sea_regions = df_enzymes[('Ocean_and_sea_regions')].unique()
enzymes = df_enzymes[('Enzyme')].unique()
enzyme_names = enzymes.tolist()

#Define color palette
color = all_palettes['Colorblind'][len(enzymes)]


# Make figure
p = figure(x_range=ocean_and_sea_regions,
           width=800,
           height=500,
           x_axis_label='Ocean and sea regions',
           y_axis_label='Total transcripts per million sequences',     
           tools="")

# Add stacked bars
p.vbar_stack(enzymes, x='Ocean_and_sea_regions', 
            width=0.9, source = source, fill_color=color,
             legend_label=enzyme_names)


#Styling
p.xaxis.axis_label_text_font_size = "16pt"
p.yaxis.axis_label_text_font_size = "16pt"
p.xaxis.major_label_text_font_size = "12pt"
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
export_png(p, filename=f'{homedir}/figures/genetics/metatranscriptomics/stab8_curson2018_bact_oceanregs.png')

# %%
