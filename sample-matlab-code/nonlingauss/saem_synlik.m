function [THETAsaem,simsummaries_final] = saem_synlik(model_param,parmask,parbase,yobs,saem_numit,warmup,numsim)
%this is the main SAEM-SL algorithm

global means_all_old_external cov_all_old_external startnegloglik simsummaries


bigtheta = model_param{1};
problem = model_param{2};
time    = model_param{3};
numdepvars = model_param{4};
vrbl    = model_param{5};

theta = param_mask(bigtheta,parmask);
fprintf('\n SAEM parameters starting values are:');
fprintf('\n %d',theta);
fprintf('\n')


THETAsaem = zeros(saem_numit+1,length(theta));
THETAsaem(1,:) = theta;

nobs = length(yobs);

[~,sobs] = feval([problem, '_summaries'],[],yobs);  % data summaries

xhat = feval([problem, '_statemodel'],bigtheta,nobs,numsim);  % X simulated
yobssim = feval([problem, '_errormodel'],bigtheta,xhat,nobs,numsim); % Y simulated
[summariesx,summariesy] = feval([problem, '_summaries'],xhat,yobssim);  % S(X) and S(Y)
summaries_all = [summariesx;summariesy];
dimx = size(summariesx,1);
% means_all and cov_all below will actually be unused. We only need to
% extract their dimensions
means_all = [mean(summariesx,2);mean(summariesy,2)];
cov_all = cov(summaries_all');

means_all_old_external = zeros(length(means_all),1); % starting means
cov_all_old_external = 1e-12*eye(size(cov_all,1));   % starting covariance

for saem_iter = 2:saem_numit+1
   fprintf('\nSAEM iteration %d',saem_iter) 
    
   % SAEM coefficients
   if saem_iter < warmup
         alpha_sequence = 1;
   else
         alpha_sequence = 1/(saem_iter-warmup+1)^0.6;
   end
   
   theta_start = THETAsaem(saem_iter-1,:);
   means_all_old = means_all_old_external;
   cov_all_old = cov_all_old_external;
   startnegloglik = inf;  % reset at each SAEM iteration
   
   % sample (summaries for the) latent state
   means_x = means_all_old(1:dimx);
   means_y = means_all_old(dimx+1:end);
   cov_x =   cov_all_old(1:dimx,1:dimx);
   cov_xy = cov_all_old(1:dimx,dimx+1:end);
   cov_yx = cov_all_old(dimx+1:end,1:dimx);
   cov_y =   cov_all_old(dimx+1:end,dimx+1:end);
  %compute conditional moments for a gaussian X|sobs
  % meansx_cond = means_x + cov_xy*inv(cov_y)*(sobs - means_y);
   meansx_cond = means_x + cov_xy*(cov_y\(sobs - means_y)); % same as the line above but more accurate
  % covx_cond = cov_x-cov_xy*inv(cov_y)*cov_yx;
   covx_cond = cov_x-cov_xy*(cov_y\cov_yx); % same as the line above but more accurate

   % we wish to generate a conditional sample.
   % attempt to generate from a multivariate Gaussian using a cholesky
   % decomposition
   [Lfactor, notposdef] = chol(covx_cond,'lower');
   if notposdef==0  % good! covx_cond is positive definite!
      sampledsumx_cond = meansx_cond + Lfactor*randn(dimx,1);
   else
       % the attempt above has failed.
       % find the nearest simmetric positive definite matrix https://se.mathworks.com/matlabcentral/fileexchange/42885-nearestspd 
       covx_cond = nearestSPD(covx_cond);
       sampledsumx_cond = mvnrnd(meansx_cond,covx_cond); % conditional sample
   end 
   
   summary_complete = [sampledsumx_cond ;sobs];

   myoptsetup = optimset('maxiter',40,'MaxFunEvals',1000);
   theta_opt = fminsearch(@(theta) negsynlik(theta,problem,parmask,parbase,nobs,numsim,summary_complete,means_all_old,cov_all_old,alpha_sequence),theta_start,myoptsetup);
   
   THETAsaem(saem_iter,:) = theta_opt;
end

%THETAmatrix_temp = [THETAmatrix_temp;THETAmatrix_saem];
save('THETAsaem','THETAsaem')

simsummaries_final = simsummaries;
save('simsummaries_final','simsummaries_final')
    
% THETAmatrix = THETAmatrix_temp;
% save('THETAmatrix','THETAmatrix')

fprintf('\n')

end


function out = negsynlik(theta,problem,parmask,parbase,nobs,numsim,summary_complete,means_all_old,cov_all_old,alpha_sequence)
   %  persistent means_all_old cov_all_old 
     global means_all_old_external cov_all_old_external startnegloglik simsummaries

      bigtheta = param_unmask(theta,parmask,parbase);
      xhat = feval([problem, '_statemodel'],bigtheta,nobs,numsim);
      yobssim = feval([problem, '_errormodel'],bigtheta,xhat,nobs,numsim);
      [summariesx,summariesy] = feval([problem, '_summaries'],xhat,yobssim);
      simsummaries = [summariesx;summariesy]; % ns x numsim matrix; ns is the total dimension of the summaries. Initial rows contain sumaries for X; then remaining rows contain summaries for Y
     % [cov_all,means_all] = robustcov(simsummaries','Method','olivehawkins');
     %  means_all = means_all';
      means_all = mean(simsummaries,2);
      cov_all = cov(simsummaries');
      
      out = neglogmvnpdf(summary_complete',means_all',cov_all);
      
      if out < startnegloglik  
          % the SAEM updating equations
          means_all_new = means_all_old + alpha_sequence*(means_all - means_all_old);
          cov_all_new   = cov_all_old   + alpha_sequence*(cov_all   - cov_all_old);
          means_all_old_external = means_all_new;  
          cov_all_old_external = cov_all_new;
          startnegloglik = out;
      end
end

