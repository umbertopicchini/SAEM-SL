# SAEM-SL
The sample-matlab-code folder contains two MATLAB implementations of the SAEM-SL methodology by Umberto Picchini, described in
"Likelihood-free stochastic approximation EM for inference in complex models", arXiv:1609.03508. 

Please refer also to the README files in each of the subfolders "nonlingauss" and "theophylline" for further information.

The software contained in the two subfolders calls the following utilities:
- neglogmvnpdf                        utility to compute the negative log pdf of a multivariate Gaussian 
- nearestSPD                          utility to compute N. Higham's "nearest semi-positive definite matrix"
- param_mask                          utility: extract the free parameters from the full vector of free and fixed parameters
- param_unmask                        utility: reconstruct the full vector of structural model parameters from the set of free parameters 
- hline utility that adds an horizontal line to plots 
