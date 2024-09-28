# CC1p0π README

The DataRelease.root file contains the data release for the manuscript "Measurement of nuclear effects in neutrino-argon interactions using generalized kinematic imbalance variables with the MicroBooNE detector"
Namely, it includes the data cross section results with the total uncertainties (denoted as eg. "TotalUnc_DeltaPn"), the covariance matrices (denoted as eg. "Cov_DeltaPn"), and the additional smearing matrix Ac (denoted as eg. "Ac_DeltaPn") by which theory predictions need to be multiplied before they are compared to the reported results. The Ac matrices need to first be applied on independent theory predictions, before the division by the bin width takes place.

The double-differential results are denoted as eg. "SerialDeltaPn_DeltaAlpha3Dq". That means that the results as a function of "DeltaPn" in slices of "DeltaAlpha3Dq" are reported in these histograms.

The cross sections are reported in units of 10^{-38} cm^{2}/Ar and the bin width division has already taken place.

The covariance matrices are reported in units of 10^{-76} cm^{4}/Ar^{2}. 
Each bin entry needs to be divided by the bin area (i.e. x bin width * y bin width) before a chi2 GoF metric is calculated.

The x-axis of the Ac matrix corresponds to the variable of interest in true space.
The y-axis of the Ac matrix corresponds to the variable of interest in regularized space.
The Ac matrices need to first be applied on independent theory predictions, before the division by the bin width takes place.

To preserve the correlations across the different phase space parts, the double-differential measurements are reported as a function of a global bin number.
The file BinScheme.txt contains the matching of the global bin number to the actual phase-space limits.
Independent theory predictions comparing to the double-differential results need to use the aforementioned global-bin numbering scheme, multiply by the corresponding Ac matrix, isolate the relevant slices of the phase-space, and finally divide by the bin width. 

