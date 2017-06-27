numattempts = 3; % how many estimation cycles? It is = 100 in the paper
all_theta_estimated = zeros(numattempts,4);

% for each trial we change the seed of the random number generator, hence
% for each traial we simulate a different dataset, using the same parameter
% values. This is the only difference between the trials

for trial=1:numattempts

rng(trial*100)

warmup = 50; % the number of SAEM iterations (K1) before using the (descreasing) "alpha" sequence into SAEM 
saem_numit = 80;   % the total number of SAEM iterations (K)
numsim = 200;  % number of simulations (R) to compute a synthetic likelihood

% Model Setup
problem = 'theophylline';                 % a string identifying the given problem at hand. 
numdepvars = 1;                           % the number of modellized "state variables"
sampletime = [1:1:30];                    % observational times
time = sampletime';  % MUST BE A COLUMN VECTOR   
stepsize = 0.05;                          % stepsize for the numerical integration of the SDE (h in the paper)
owntime = [0:stepsize:sampletime(end)];   % the fine time discretization for the numerical solution of the SDE

fprintf('\n\nCOMPUTING...\n\n');



%:::::::::::: HERE GO TRUE PARAMETER VALUES FOR DATA GENERATION :::::::::::::::
X0 = 8;
log_Ke = -3;  % Ke = 0.0498
log_Ka = 0.4;   % Ka = 1.49
log_Cl = -3.22;  % Cl = 0.04
sigma = 0.1;  % log_sigma = -2.3
sigmaepsilon = 0.32; % log_sigmaepsilon = -1.14
log_sigma = log(sigma);
log_sigmaepsilon = log(sigmaepsilon);


vrbl = [1:numdepvars];
% Now replicate the values in vrbl depending on the time value 
vrbl = repmat(vrbl,1,length(unique(time)));


%::::::::: GENERATE DATA ::::::::::::::::::::::::::::::::::::::::

bigtheta_true = [X0,log_Ke,log_Ka,log_Cl,log_sigma,log_sigmaepsilon];   % store here all parameters needed for SDE simulation

% create 'observations' by adding measurement noise to the noise-free
% trajectory
xhat_true = feval([problem, '_statemodel'],bigtheta_true,[],1,sampletime,owntime);
yobs = feval([problem, '_errormodel'],bigtheta_true,xhat_true,length(sampletime));

% plot(sampletime,yobs,'o-')
% return
% %::::::::: END OF DATA GENERATION ::::::::::::::::::::::::::::::::::::

%                 X0    log_Ke   log_Ka   log_Cl      log_sigma    log_sigmaepsilon
bigtheta_start = [8,      -1.9,      0.4,    -2,              -2,          -0.69 ];
parmask        = [0 ,      1,        0,       1,               1,            1  ];
parbase = bigtheta_start;


   model_param = {bigtheta_start,problem,owntime,time,numdepvars,vrbl};  
   covariates = {sampletime',owntime'};

   THETAmatrix_saem = saem_synlik(model_param,parmask,parbase,yobs,covariates,saem_numit,warmup,numsim);

   all_theta_estimated(trial,:) = THETAmatrix_saem(end,:);
 
    if trial==1
        figure
    end
   if isreal(THETAmatrix_saem) && any(any(isnan(THETAmatrix_saem)))==0
   % exponentiate so we transform log-parameters to their natural scale
   subplot(2,2,1)
   plot(exp(THETAmatrix_saem(:,1)))
   xlabel('Ke')
   hline(exp(log_Ke))  % line showing the true parameter value
   hold on
   subplot(2,2,2)
   plot(exp(THETAmatrix_saem(:,2)))
   xlabel('Cl')
   hline(exp(log_Cl))  % line showing the true parameter value
   hold on
   subplot(2,2,3)
   plot(exp(THETAmatrix_saem(:,3)))
   xlabel('\sigma')
   hline(exp(log_sigma))  % line showing the true parameter value
   hold on
   subplot(2,2,4)
   plot(exp(THETAmatrix_saem(:,4)))
   xlabel('\sigma_{\epsilon}')
   hline(exp(log_sigmaepsilon))  % line showing the true parameter value
   hold on
   end
   saveas(gcf,'saem_trajectories')
end
save('all_theta_estimated','all_theta_estimated')
hold off


