function yobssim = theophylline_errormodel(bigtheta,xhat,nobs)
% produce Y, i.e. add measurement noise to X

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

