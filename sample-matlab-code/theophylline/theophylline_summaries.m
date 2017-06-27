function [summariesx,summariesy] = theophylline_summaries(xhat,xhat_big,yobssim,covariates)
% compute summary statistics for both X and Y across all simulations
% (vectorised code)

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

elseif isempty(xhat) && isempty(xhat_big) 
    summariesx = [];
    summariesy = [
                  median(yobssim,1);...
                  mad(yobssim,0,1);...
                  (yobssim(end,:)-yobssim(1,:))./(sampletime(end)-sampletime(1))];
end



