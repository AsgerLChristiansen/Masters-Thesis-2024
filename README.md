# Masters-Thesis-2024
Repository with code used in the Master's Thesis of Asger L. Christiansen at the Master's Degree in Cognitive Science, Aarhus University, 2024

## Dependencies
Software dependencies are detailed in the report. Code assumes pacman is installed, as pacman is used to load / install any other packages.

# Repo Structure

## Root directory
Contains a copy of the final report, as well as the Data and Pipeline folders.

### Data Folder
The Data folder is intended to contain the raw fNIRS data, and exists mainly for illustrative purposes, as the data is not present in this repository. Any and all files created after loading the raw data (including time series of raw haemoglobin concentration, WTC analysis output, et cetera) are stored in subdirectories of the Pipeline folder.

### Pipeline Folder
The Pipeline folder includes several numbered scripts and notebooks, some Python-based, some R-based, as well as the various folders those scripts use. Empty folders are included, again, for illustrative purposes. The first script, "1. MNE Preprocessing.ipynb", will naturally not run without the raw fNIRS data, and subsequent scripts and notebooks require that script to have run first.
As a result, do not expect to be able to run the code. The present README will provide a brief overview of what each script does, in numbered order.

Number 1 is a Jupyter Notebook that takes the raw fNIRS data and extracts time series of oxy- and deoxyhaemoglobin concentration from them, and organizes them along pairs, channels and visits in the Pipeline/HaemoglobinTimeSeries folder. This prepares the time series for WTC analysis.

Number 2 creates the 22 surrogate pairs used in the study.

Number 3 and 4 perform wavelet transform coherence (WTC) analysis on both real and surrogate pairs, for every channel and visit, and stores the many files that this method outputs.

Number 5 computes average coherence values during each experimental condition (excluding the cone of influence) from all these various files, and organizes it all into a single .csv file of 11616 data points for subsequent stan analysis.

Number 6 makes a few initial observations about the data to help with prior choice, pads the dataset for missing data, and otherwise prepares it for analysis.

Number 7 runs parameter recovery for the two stan models in the Pipeline/Stan Code folder. *NOTE*: This script should be able to be run on its own, if the reader wishes to replicate my parameter recovery results on simulated data. The same seed is set in that script as was used on my machine, so results SHOULD be identical. However, it should be noted that, at least in Model 2's case, the parameter recovery in the report was made with a bug in the simulation script which has since been fixed.
If you rerun my parameter recovery, beware that it it very time-intensive! You should find that Model 1 can recover parameters reasonably, while Model 2 fails completely. Results of Parameter Recovery will be saved in the Parameter Recovery folder.

Number 8 runs the actual stan analysis of Model 1 (Since Model 2's parameter recovery was very bad), and saves results in the appropriately named folder.

Finally, Number 9 makes most of the plots already present in the Plots folder.

Folders /src and /temp contain, respectively scripts with helper functions and temporary .csv files for use in various stages of the analysis. The folder /DAGs should contain the Directed Acyclical Graphs used in the report.
