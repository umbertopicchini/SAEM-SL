function [summariesx,summariesy] = theophylline_summaries(xhat,xhat_big,yobssim,covariates)
% compute summary statistics for both X and Y across all simulations (vectorised code)
% xhat: a nobs x (dx*numsim) matrix, where nobs = #observations, dx = size of the unobserved system X (here dy=1), numsim=number of simulations from the model.
% xhat_big: a N x (dx*numsim) matrix, where N = #imputed time points (=length(owntime) in the run file), dx = size of the unobserved system X (here dy=1), numsim=number of simulations from the model.
% yobssim: a nobs x (dy*numsim) matrix, where nobs = #observations, dy = size of the observed system Y (here dy=1), numsim=number of simulations from the model.
% covariates: a structure, with first component the observational times (sampletime), and second component the fine grid of imputed times (owntime).
% summariesx = a dsx * numsim matrix, where dsx is the number of summary statistics for X.
% summariesy = a dsy * numsim matrix, where dsy is the number of summary statistics for Y.

sampletime = covariates{1};
owntime = covariates{2};
dt=owntime(2)-owntime(1);


if ~isempty(xhat) && ~isempty(xhat_big)
   % remove simulations leading to negative trajectories
   indexcomplex = find(any(xhat_big<0));
   numsim = size(xhat,2);
   if ~isempty(indexcomplex)
     ok_include = ~ismember([1:numsim],indexcomplex);
     xhat_big = xhat_big(:,ok_include);
     xhat = interp1(owntime',xhat_big,sampletime');
     yobssim = yobssim(:,ok_include);
   end 
    

   summariesx = [median(xhat_big,1);...
                mad(xhat_big,0,1);...
                sqrt(sum( (xhat_big(2:end,:)-xhat_big(1:end-1,:)).^2,1) ./ (sum(sqrt(xhat_big(1:end-1,:))*dt,1)));...
                sqrt(sum((yobssim-xhat).^2,1)/length(sampletime))];
  summariesy = [median(yobssim,1);...
                mad(yobssim,0,1);...
                (yobssim(end,:)-yobssim(1,:))./(sampletime(end)-sampletime(1))];

elseif isempty(xhat) && isempty(xhat_big) % only used when computing summaries for observed data
    summariesx = [];
    summariesy = [
                  median(yobssim,1);...
                  mad(yobssim,0,1);...
                  (yobssim(end,:)-yobssim(1,:))./(sampletime(end)-sampletime(1))];
end



