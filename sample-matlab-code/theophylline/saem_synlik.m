function THETAsaem = saem_synlik(model_param,parmask,parbase,yobs,covariates,saem_numit,warmup,numsim)
% this is the main SAEM-SL algorithm (notice this version is not very general. It is specialised for the specific example "theophylline").
% model_param: a structure containing quantities defined in the run file (bigtheta,problem,owntime,sampletime,numdepvars,vrbl). See theophylline_run.m.
% parmask: vector containing 1's for free parameters and 0 otherwise, see See theophylline_run.m.
% parbase: contains starting values for the parameters (especially useful to extract the values of fixed parameters = not estimated parameters). See See theophylline_run.m.
% yobs: vector of observations from the Y component of the model.
% covariates: ("covariates" is redundant and should be eliminated) a structure, with first component the observational times (sampletime), and second component the fine grid of imputed times (owntime).
% saem_numit: total number of SAEM iterations (=K in the paper).
% warmup: the number of SAEM iterations (K1) before using the (descreasing) "alpha" sequence into SAEM
% numsim: number of simulations (R in the paper) to compute a single evaluation of a synthetic likelihood.
% THETAsaem = a (sam_numit * sum(parmask)) matrix containing the evolution of the estimated parameters.


global means_all_old_external cov_all_old_external startnegloglik


bigtheta = model_param{1};
problem = model_param{2};
owntime    = model_param{3};
sampletime    = model_param{4};
numdepvars = model_param{5};
vrbl    = model_param{6};

theta = param_mask(bigtheta,parmask);
fprintf('\n SAEM parameters starting values are:');
fprintf('\n %d',theta);
fprintf('\n')


THETAsaem = zeros(saem_numit+1,length(theta));
THETAsaem(1,:) = theta;

nobs = length(yobs);

[~,sobs] = feval([problem, '_summaries'],[],[],yobs,covariates);

[xhat,xhat_big] = feval([problem, '_statemodel'],bigtheta,nobs,numsim,sampletime,owntime);
yobssim = feval([problem, '_errormodel'],bigtheta,xhat,nobs);
[summariesx,summariesy] = feval([problem, '_summaries'],xhat,xhat_big,yobssim,covariates);
summaries_all = [summariesx;summariesy];
dimx = size(summariesx,1);
% means_all and cov_all below will actually be unused. We only need to
% extract their dimensions
means_all = [mean(summariesx,2);mean(summariesy,2)];
cov_all = cov(summaries_all');

%initialize mean and covariance
means_all_old_external = zeros(length(means_all),1);
cov_all_old_external = 1e-12*eye(size(cov_all,1));

for saem_iter = 2:saem_numit+1
   fprintf('\nSAEM iteration %d',saem_iter) 
   
   if saem_iter < warmup
         alpha_sequence = 1;
   else
         alpha_sequence = 1/(saem_iter-warmup+1);
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
  % covx_cond = (covx_cond + covx_cond.') / 2;  % same as covx_cond. Just a trick to avoid the mvnrnd error "Error using mvnrnd SIGMA must be a symmetric positive semi-definite matrix."
  % sampledsumx_cond = mvnrnd(meansx_cond,covx_cond); % conditional sample
   
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

   myoptsetup = optimset('maxiter',30,'MaxFunEvals',10000);
   theta_opt = fminsearch(@(theta) negsynlik(theta,problem,parmask,parbase,nobs,sampletime,owntime,covariates,numsim,summary_complete,means_all_old,cov_all_old,alpha_sequence),theta_start,myoptsetup);

   THETAsaem(saem_iter,:) = theta_opt;
   THETAsaem_temp = THETAsaem(1:saem_iter,:);
   save('THETAsaem_temp','THETAsaem_temp')  
end

%THETAmatrix_temp = [THETAmatrix_temp;THETAmatrix_saem];
save('THETAsaem','THETAsaem')

    
% THETAmatrix = THETAmatrix_temp;
% save('THETAmatrix','THETAmatrix')

fprintf('\n')

end


function out = negsynlik(theta,problem,parmask,parbase,nobs,sampletime,owntime,covariates,numsim,summary_complete,means_all_old,cov_all_old,alpha_sequence)
   %  persistent means_all_old cov_all_old 
     global means_all_old_external cov_all_old_external startnegloglik
     
      bigtheta = param_unmask(theta,parmask,parbase);
      [xhat,xhat_big] = feval([problem, '_statemodel'],bigtheta,nobs,numsim,sampletime,owntime);
      yobssim = feval([problem, '_errormodel'],bigtheta,xhat,nobs);
      [summariesx,summariesy] = feval([problem, '_summaries'],xhat,xhat_big,yobssim,covariates);
      simsummaries = [summariesx;summariesy]; % ns x numsim matrix; ns is the total dimension of the summaries. Initial rows contain sumaries for X; then remaining rows contain summaries for Y
      [cov_all,means_all] = robustcov(simsummaries','Method','olivehawkins');
   %   [cov_all,means_all] = robustcov(simsummaries');
   
      means_all = means_all';
    %  means_all = mean(simsummaries,2);
    %  means_all = means_all';
    %  cov_all = cov(simsummaries');
      
      
      out = neglogmvnpdf(summary_complete',means_all',cov_all);
      
      if out < startnegloglik
          means_all_new = means_all_old + alpha_sequence*(means_all - means_all_old);
          cov_all_new   = cov_all_old   + alpha_sequence*(cov_all   - cov_all_old);
          means_all_old_external = means_all_new;
          cov_all_old_external = cov_all_new;
          startnegloglik = out;
      end
end

