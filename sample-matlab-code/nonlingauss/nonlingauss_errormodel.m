function yobssim = nonlingauss_errormodel(bigtheta,xhat,nobs,numsim)

% produce Y, i.e. add measurement noise to X

X0 = bigtheta(1);
log_sigmax = bigtheta(2);
log_sigmay = bigtheta(3);

sigmay = exp(log_sigmay);

yobssim = xhat + sigmay*randn(nobs,numsim);


end

