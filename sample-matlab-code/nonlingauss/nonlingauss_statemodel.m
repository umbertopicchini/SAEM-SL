function xhat = nonlingauss_statemodel(bigtheta,nobs,numsim)
% produce nobs x numsim simulations of the X coordinate of the model
% nobs = number of sampling times
% numsim = number of simulations (R in the paper)

X0 = bigtheta(1);
log_sigmax = bigtheta(2);
log_sigmay = bigtheta(3);

sigmax = exp(log_sigmax);

xhat = zeros(nobs,numsim);

xhat(1,:) = 2*sin(exp(X0))*ones(1,numsim) + sigmax*randn(1,numsim);
xhat_pre = xhat(1,:);

for ii=2:nobs
    xhat(ii,:) = 2*sin(exp(xhat_pre)) + sigmax*randn(1,numsim);
    xhat_pre = xhat(ii,:);
end

end

