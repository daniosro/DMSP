## `data`

This folder contains all the raw and processed data from the project. It contains the following subdirectories:

 * **`raw`**: Includes all the raw data, subdivided as:

  1. EA: Sulfur isotope data from the EA-IRMS for the enzymatic degradation of DMSP by DMSP lyases and demethylases.
  2. HPLC: DMSP composition from HPLC for the enzymatic degradation of DMSP by DMSP lyases and demethylases.
  3. enz_deg: DMSP composition from HPLC by the eukaryotic DMSP lyase Alma1, in experiments targeted to assess loss of enzyme activity.
  4. genetics: Metagenomics and metatranscriptomics data from the Tara Oceans expedition, as published by Curson *et al.* (2018) and Landa *et al.* (2019).

 * **`processed`**: Includes the processed data, subdivided as:

  1. EA: Sulfur isotope data from the EA-IRMS, corrected by blank, linearity and standards.
  2. HPLC: DMSP composition from HPLC and fraction of DMSP remaining at each timepoint for all the enzymes. 
  3. genetics: Tidy dataframes from metagenomics and metatranscriptomics data from the Tara Oceans expedition, as published by Curson *et al.* (2018) and Landa *et al.* (2019).
  4. enzymes: Corrected data of degradation and sulfur isotopic composition of DMSP for each enzyme degradation experiment.
  
 * **`modelling`**: Includes the modelling generated data used to produce the Supp. Figure 4 (changes in the fractionation factor of DddP with different assumptions of enzyme degradation rate, V<sub>max</sub> and V<sub>max</sub>/K<sub>M</sub> and the Figure 4 (variations in the &delta;<sup>34</sup>S of seawater DMSP affected by demethylation, bacterial cleavage and eukaryotic cleavage) of the paper.
