### On-/Off-Axis Data Release
#### (Version 1.0.1, dated 2024/08/12)

This tar archive contains the data release for ‘First measurement of muon neutrino charged-current interactions on hydrocarbon without pions in the final state using multiple detectors with correlated energy spectra at T2K’. It contains the cross-section data points and supporting information in ROOT and text format, which are detailed below:

+ `onoffaxis_xsec_data.root`
This ROOT file contains the extracted cross section and the nominal MC prediction as TH1D histograms for both the flattened 1D array of bins and in the angle binning for the analysis. The ROOT file also contains both the covariance and inverted covariance matrix for the result stored as TH2D histograms. The angle bin numbering and the corresponding bin edges are detailed at the end of the README.

+ `flux_analysis.root`
This ROOT file contains the nominal and post-fit flux histograms for ND280 and INGRID. Two different binnings are included: a fine binned histogram (220 bins) and a coarse binned histogram (20 bins). The coarse binned histogram corresponds to the flux parameters detailed in the paper (and bin edges listed in the appendix).

+ `xsec_data_mc.csv`
The extracted cross-section data points and the nominal MC prediction for each bin is stored as a comma-separated value (CSV) file with header row.

+ `cov_matrix.csv` and `inv_matrix.csv`
The covariance matrix and the inverted covariance matrix are both stored as CSV files with each row stored as a single line and columns separated by commas (there is no header row). Matrix element (0,0) corresponds to the first number in the file.

+ `nd280_analysis_binning.csv` and `ingrid_analysis_binning.csv`
The analysis bin edges are included as CSV files. The columns are labeled with a header row and denote the linear bin index and the lower and upper bin edge for the angle and momentum bins. The units are in cos(angle) for the angle bins and in MeV/c for the momentum bins.

+ `calc_chisq.cxx`
This is an example ROOT script to calculate the chi-square between the data and the nominal MC prediction using the ROOT file in the data release. To run, open ROOT and load the script (`.L calc_chisq.cxx`) and execute the function `calc_chisq("/path/to/file.root")`.

+ `calc_chisq.py`
This is an example Python script to calculate the chi-square between the data and the nominal MC prediction using the text/CSV files in the data release. The code requires NumPy as an external dependency, but otherwise uses built-in modules. To run, execute using a Python3 interpreter and give the file paths to the data/MC text file and the inverse covariance text file as the first and second arguments respectively -- e.g. `python3 calc_chisq.py /path/to/xsec_data_mc.csv /path/to/inv_matrix.csv`

+ ND280 angle bin numbering
    - 0: `-1.0 < cos(#theta) < 0.20`
    - 1: `0.20 < cos(#theta) < 0.60`
    - 2: `0.60 < cos(#theta) < 0.70`
    - 3: `0.70 < cos(#theta) < 0.80`
    - 4: `0.80 < cos(#theta) < 0.85`
    - 5: `0.85 < cos(#theta) < 0.90`
    - 6: `0.90 < cos(#theta) < 0.94`
    - 7: `0.94 < cos(#theta) < 0.98`
    - 8: `0.98 < cos(#theta) < 1.00`

+ INGRID angle bin numbering
    - 0: `0.50 < cos(#theta) < 0.82`
    - 1: `0.82 < cos(#theta) < 0.94`
    - 2: `0.94 < cos(#theta) < 1.00`
    
### Changelog

#### v1.0.1
Fix transcription error in INGRID momentum binning. The lowest momentum bin edge is at 350 MeV/c, not 300 MeV/c.
