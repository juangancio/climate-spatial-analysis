# climate-spatial-analysis

Python scripts for spatial analysis using ordinal patterns in climate data.
Code is preprared for analysing two regions (El Nino 3.4 and Gulf Stream) from two data sets (NOAA OI v2 and ERA5 global reanalysis). But can be easily generalised for any input data.

spe_analysis.py: Loads the data, calculates the anomaly, extracts the Spatial Permutation Entropy (SPE) from the data, in the horizontal (W-E) direction and in the vertical (N-S) direction, saves these signals, and plots them.

spe_utils.py: Utilities for SPE analysis.

Study:trend.py: Script for significance test of signals' trends.

make_figs.m: Matlab scrip for figures.

