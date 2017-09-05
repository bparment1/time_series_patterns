# time_series_patterns
This repository contains general code to analyze time series. This is as part of Research Support at SESYNC.
The goal is to provide a set of functions related analysis of time series, trends and variability and cycles.

This includes:
1) detection and removeal of trends:
- theil sen
- ordinary least square
- differencing
- moving averages/kernels
2) identification and removal of periodic cycles
- fast fourier transform spectrum
- fast fourier transform periodogram
- multitaper methods
- windowed fourier transform
- wavelets
- PCA on lag analyses
3) Variability
- PCA S-mode
- PCA T-mode
- Space-time pca/eof on spatio-temporal data
4) General tools
- detection of autocorrelation
- significance tests for various methods

