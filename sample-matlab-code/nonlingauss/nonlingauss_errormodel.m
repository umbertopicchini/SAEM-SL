function yobssim = nonlingauss_errormodel(bigtheta,xhat,nobs,numsim)

% simulates Y, i.e. adds measurement noise to X.
%
% bigtheta = vector of the model parameters (both free and fixed parameters).
% xhat = a nobs x numsim matrix of of simulations from the latent state X
% nobs = number of sampling times
% numsim = number of simulations (=R in the paper).
% yobssim = a nobs x numsim matrix of of simulations from the observable process Y


X0 = bigtheta(1);
log_sigmax = bigtheta(2);
log_sigmay = bigtheta(3);

sigmay = exp(log_sigmay);

yobssim = xhat + sigmay*randn(nobs,numsim);


end

