## `data`

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