function [yobs,xhat] = theophylline_modelsimulate(bigtheta,time)

% Exact simulation of the solution of the SDE model for the Theophylline example.
%
% BIGTHETA is the complete structural parameter vector.
% OWNTIME is the time-grid for the numerical discretization.
% NUMDEPVARS is the size of the sde model (=1 in this case).
% XHAT is the solution of the SDE model at time-grid OWNTIME.

% This file is part of the "abc-sde" program. Copyright (C) 2013 Umberto Picchini
% https://sourceforge.net/projects/abc-sde/
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Parameters
xzero     = bigtheta(1);
log_Ke    = bigtheta(2);
log_Ka    = bigtheta(3);
log_Cl    = bigtheta(4);
log_sigma = bigtheta(5);
log_sigmaepsilon = bigtheta(6);


Ke = exp(log_Ke);
Ka = exp(log_Ka);
Cl = exp(log_Cl);
sigma = exp(log_sigma);
sigmaepsilon = exp(log_sigmaepsilon);

Dose = 4;  % the "dose"
% generate a trajectory using the SDE EXACT solution
nobs = length(time);
xhat = zeros(nobs,1);
W = [0,cumsum(sqrt(time(2)-time(1))*randn(1,nobs-1))];
ito_sum_integrand = arrayfun(@(index) exp(Ke*time(index-1))*(W(index)-W(index-1)),2:nobs);
ito_sum_integrand = [0, ito_sum_integrand];
ito_sum = cumsum(ito_sum_integrand);
for i = 1:nobs
    xhat(i) = xzero*exp(-Ke*time(i)) + Dose*Ka*Ke/(Cl*(Ke-Ka)) * (exp(-Ka*time(i))-exp(-Ke*time(i))) + sigma*exp(-Ke*time(i)).* ito_sum(i);
end

yobs = xhat + sigmaepsilon*randn(nobs,1);

end

