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
# Make stacked bar plot for marine pelagic biomes following 
# the Longhurst code: https://en.wikipedia.org/wiki/Longhurst_code. 

# Get averages for the Percentage of genome equivalents for each marine pelagic biome and enzyme
df_stab3_landaetal_2019_enzymes_averages2 = df_stab3_landaetal_2019_enzymes.groupby(['Marine_pelagic_biomes _Longhurst_2007', 'Gene']).Percentage_of_genome_equivalents.agg(['mean','std']).reset_index()
#Pivot table to organize the means by enzyme and marine pelagic biome
df_stab3_landaetal_2019_enzymes_averages2 = df_stab3_landaetal_2019_enzymes_averages2.pivot_table(index=['Marine_pelagic_biomes _Longhurst_2007'], columns='Gene', values='mean').reset_index()
df_stab3_landaetal_2019_enzymes_averages2.head()
# %%
#Make stacked plot for enzymes with depth zones

#Define series and lists of unique values
enzymes = df_stab3_landaetal_2019_enzymes_averages[('Gene')].unique()
enzyme_names = enzymes.tolist()
marine_pelagic_biomes = df_stab3_landaetal_2019[('Marine_pelagic_biomes _Longhurst_2007')].unique()

#Make stacked plot for marine pelagic biomes with enzymes
#Define data source
source = ColumnDataSource(data=df_stab3_landaetal_2019_enzymes_averages2)
color = all_palettes['Colorblind'][len(enzymes)]

# Make figure
p = figure(x_range=marine_pelagic_biomes,
           width=800,
           height=500,
           x_axis_label='Marine pelagic biomes',
           y_axis_label='Percentage of genome equivalents',     
           tools="")

# Add stacked bars
p.vbar_stack(enzymes, x='Marine_pelagic_biomes _Longhurst_2007', 
            width=0.9, source = source, fill_color=color,
             legend_label=enzyme_names)

#Styling
p.xaxis.axis_label_text_font_size = "16pt"
p.yaxis.axis_label_text_font_size = "16pt"
p.xaxis.major_label_text_font_size = "7pt"
p.yaxis.major_label_text_font_size = "14pt"


bokeh.io.show(p)
#Export figure
export_png(p, filename=f'{homedir}/figures/genetics/metagenomics/stab3_landa2019_bact_biomes.png')
# %%
