function [xhat,xhat_big] = theophylline_statemodel(bigtheta,nobs,numsim,sampletime,owntime)

% Produces am N x numsim matrix xhat_big of simulations of the X coordinate
% of the model, where N=length(owntime).
% Also produces an nobs x numsim matrix xhat of simulations of the X coordinate
% of the model, where nobs = length(sampletime). This are linear
% interpolated values from the xhat_big output.
% nobs = number of sampling times
% numsim = number of simulations (R in the paper)

% Parameters
xzero     = bigtheta(1);
log_Ke    = bigtheta(2);
log_Ka    = bigtheta(3);
log_Cl    = bigtheta(4);
log_sigma = bigtheta(5);

Ke = exp(log_Ke);
Ka = exp(log_Ka);
Cl = exp(log_Cl);
sigma = exp(log_sigma);

Dose = 4; % the adminstered drug dose

N= length(owntime);
xhat_big = zeros(N,numsim);

xhat_big(1,:) = xzero;
xhat_pre = xhat_big(1,:);

for j=2:N

    h = owntime(j)- owntime(j-1);
    t = owntime(j-1);
  
    Winc = sqrt(h)*randn(1,numsim); % the Wiener increment(s) dWj;
    driftX = (Dose * Ka*Ke)/Cl * exp(-Ka*t) -Ke * xhat_pre ;          % the Ito SDE drift
    diffusionX = sigma*sqrt(abs(xhat_pre));      % the Ito SDE diffusion

    xhat_next = xhat_pre + driftX * h + diffusionX .* Winc;  
    xhat_big(j,:) = xhat_next;
    xhat_pre = xhat_next;
end



xhat = interp1(owntime',xhat_big,sampletime');

end


