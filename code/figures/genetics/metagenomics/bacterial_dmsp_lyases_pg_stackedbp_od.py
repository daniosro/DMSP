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

# Load data
df_stab3_landaetal_2019 = pd.read_csv(f'{homedir}/data/processed/genetics/stab3_landaetal_2019_tidy.csv')
df_stab3_landaetal_2019.head()

# %%
#Data treatment to make the plot for the percentage of genome equivalents for each enzyme

# Get the total number of genome equivalents
#Filter by total
df_reca = df_stab3_landaetal_2019[(df_stab3_landaetal_2019.Gene == 'recA_hits')]

#Filter out the recA hits
df_stab3_landaetal_2019_enzymes = df_stab3_landaetal_2019[df_stab3_landaetal_2019.Gene != 'recA_hits'] 

# Create dataframe for the mean percentage of genome equivalents for each enzyme
df_stab3_landaetal_2019_enzymes_averages = df_stab3_landaetal_2019_enzymes.groupby('Gene').Percentage_of_genome_equivalents.agg(['mean','std']).reset_index()
#Show dataframe
df_stab3_landaetal_2019_enzymes_averages

# %%
# Make stacked bar plot for ocean depths.
# DCM = Depth of chlorophyll maximum, 
# MES = mesopelagic ocean, MIX = mixed layer depth, 
# SRF = surface ocean. 

# Get averages for the total number of transcripts per million sequences for each enzyme and depth zone
df_stab3_landaetal_2019_enzymes_averages = df_stab3_landaetal_2019_enzymes.groupby(['Gene', 'Depth_zone']).Percentage_of_genome_equivalents.agg(['mean','std']).reset_index()
#Pivot table to organize the means by enzyme and depth zone
df_stab3_landaetal_2019_enzymes_averages = df_stab3_landaetal_2019_enzymes_averages.pivot_table(index=['Gene'], columns='Depth_zone', values='mean').reset_index()
df_stab3_landaetal_2019_enzymes_averages
# %%
#Make stacked plot for enzymes with depth zones

#Define series and lists of unique values
enzymes = df_stab3_landaetal_2019_enzymes_averages[('Gene')].unique()
depth_zones = df_stab3_landaetal_2019[('Depth_zone')].unique()
depth_zone_names = depth_zones.tolist()

#Define data source
source = ColumnDataSource(data=df_stab3_landaetal_2019_enzymes_averages)
color = all_palettes['Colorblind'][len(depth_zones)]

# Make figure
p = figure(x_range=enzymes,
           width=800,
           height=500,
           x_axis_label='Gene',
           y_axis_label='Percentage of genomes with each gene',     
           tools="")

# Add stacked bars
p.vbar_stack(depth_zones, x='Gene', 
            width=0.9, source = source, fill_color=color,
             legend_label=depth_zone_names)

#Styling
p.xaxis.axis_label_text_font_size = "16pt"
p.yaxis.axis_label_text_font_size = "16pt"
p.xaxis.major_label_text_font_size = "14pt"
p.yaxis.major_label_text_font_size = "14pt"


bokeh.io.show(p)
#Export to png
export_png(p, filename=f'{homedir}/figures/genetics/metagenomics/stab3_landa2019_bact_depths.png')
# %%
