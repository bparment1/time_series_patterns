# time_series_patterns
This repository contains general code to analyze time series. This is generated as part of the Research Support at SESYNC.
The goal is to provide a set of functions related to the analysis of time series with emphasis on trends, variability and cycles.

This includes:
1) detection and removal of trends:
- theil sen
- ordinary least square
- differencing
- moving averages/kernels
2) identification and removal of periodic cycles
- harmonic regression
- fast fourier transform spectrum
- fast fourier transform periodogram
- multitaper methods
- windowed fourier transform/ Short Term Fourier Transform
- wavelets
- PCA on lag analyses
- filtering of frequencies using window function and FIR filters
3) Variability
- PCA S-mode
- PCA T-mode
- Space-time pca/eof on spatio-temporal data
4) General tools
- detection of autocorrelation
- significance tests for various methods
- generation of synthetic time series with trends, cycles and random components

