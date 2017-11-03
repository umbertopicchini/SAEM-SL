function yobssim = theophylline_errormodel(bigtheta,xhat,nobs)
% simulates Y, i.e. adds measurement noise to X.
% 
% bigtheta = vector of the model parameters (both free and fixed parameters).
% xhat = a nobs x numsim matrix of of simulations from the latent state X
% nobs = number of sampling times
% numsim = number of simulations (=R in the paper).
% yobssim= a nobs x numsim matrix of of simulations from the observable process Y

X0 = bigtheta(1);
log_Ke = bigtheta(2);
log_Ka = bigtheta(3);
log_Cl = bigtheta(4);
log_sigma = bigtheta(5);
log_sigmaepsilon = bigtheta(6);

sigmaepsilon = exp(log_sigmaepsilon);
numsim = size(xhat,2);
yobssim = xhat + sigmaepsilon*randn(nobs,numsim);


end

