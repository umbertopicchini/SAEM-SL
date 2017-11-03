function [summariesx,summariesy] = nonlingauss_summaries(xhat,yobssim)
% compute summary statistics for both X and Y across all simulations (vectorised code)
% xhat: a nobs x (dx*numsim) matrix, where nobs = #observations, dx = size of the unobserved system X (here dy=1), numsim=number of simulations from the model.
% xhat: a nobs x (dy*numsim) matrix, where nobs = #observations, dy = size of the observed system Y (here dy=1), numsim=number of simulations from the model.


if ~isempty(xhat)
   summariesx = [
                median(xhat,1);...
                mad(xhat,1,1);...
                prctile(xhat,10);...
                prctile(xhat,25);...
                prctile(xhat,75);...
                prctile(xhat,90)];
  summariesy = [
                median(yobssim,1);...
                mad(yobssim,1,1);...
                prctile(yobssim,10);...
                prctile(yobssim,25);...
                prctile(yobssim,75);...
                prctile(yobssim,90)];
elseif isempty(xhat)  % used only to compute summaries of the observed data
    summariesx = [];
    summariesy = [
                  median(yobssim,1);...
                  mad(yobssim,1,1);
                  prctile(yobssim,10);...
                  prctile(yobssim,25);...
                  prctile(yobssim,75);...
                  prctile(yobssim,90)];
end


end

