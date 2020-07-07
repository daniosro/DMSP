# Sulfur isotope fractionations constrain the biological cycling of dimethylsulfoniopropionate (DMSP) in the Upper Ocean

Welcome to the GitHub repository that includes the code and data for this project.

## Repository Architecture

The repository is split into three main directories. Please see each directory for information about each file.

### **`code`** 
This folder includes all the source code used for processing the data, doing the modelling, and generating all figures. It is divided in the following subdirectories:
 * **`processing`**: Includes all the code used to correct and transform the data. 
 * **`analysis`**: Includes the scripts used to make Michaelis-Menten fits, linear regressions, and other analysis of the data. 
 * **`modelling`**: Contains the scripts used to model the changes in the isotopic composition of DMSP in vitro (under a scenario of enzyme degradation) and in the environment (with multiple degradation pathways).
 * **`figures`**: Contains all the scripts used to generate the figures in the main text and the supplementary material. 
  * **`notebooks`**: Contains all the jupyter notebooks that describe the data treatment, modelling, and figure plotting. 

 ### **`templates`** 
This folder contains templates for all the analysis performed in this project.

### **`data`** 
This folder contains all the raw and processed data from the project. It contains the following subdirectories:

 * **`raw`**: Includes all the raw data, subdivided as:

  1. EA: Sulfur isotope data from the EA-IRMS for the enzymatic degradation of DMSP by DMSP lyases and demethylases.
  2. HPLC: DMSP composition from HPLC for the enzymatic degradation of DMSP by DMSP lyases and demethylases.
  3. enz_deg: DMSP composition from HPLC by the eukaryotic DMSP lyase Alma1, in experiments targeted to assess loss of enzyme activity.
  4. genetics: Metagenomics and metatranscriptomics data from the Tara Oceans expedition, as published by Curson *et al.* (2018) and Landa *et al.* (2019).

 * **`processed`**: Includes the processed data, subdivided as:

  1. EA: Sulfur isotope data from the EA-IRMS, corrected by blank, linearity and standards.
  2. genetics: Tidy dataframes from metagenomics and metatranscriptomics data from the Tara Oceans expedition, as published by Curson *et al.* (2018) and Landa *et al.* (2019).
  3. enzymes: Corrected data of degradation and sulfur isotopic composition of DMSP for each enzyme degradation experiment.
 
# License Information

All creative works (code, figures, etc) are licensed under the [Creative
Commons CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) license. This work is published from: United States.
