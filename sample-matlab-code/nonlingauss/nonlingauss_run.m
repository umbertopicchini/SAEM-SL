
rng(10)

numattempts = 30;  % how many estimation cycles? It is = 30 in the paper
warmup = 10; % the number of SAEM iterations (K1) before using the (descreasing) "alpha" sequence into SAEM 
saem_numit = 20;   % the total number of SAEM iterations (K)
numsim = 500;  % number of simulations (R) to compute a synthetic likelihood


% Model Setup
problem = 'nonlingauss';                 % a string identifying the given problem at hand. 
numdepvars = 1;                           % the number of modelled "state variables"
integrator = [];                    % keep this empty
sampletime = [1:1:50];              % observational times
time = sampletime';  % MUST BE A COLUMN VECTOR              
owntime = time;  % NOT NECESSARY, WILL NOT BE USED


fprintf('\n\nCOMPUTING...\n\n');



%:::::::::::: HERE GO TRUE PARAMETER VALUES FOR DATA GENERATION :::::::::::::::
X0 = 0;  % initial state
sigmax = sqrt(5);  % = 2.23 -->  % log_sigmax = 0.805
sigmay = sqrt(5); %  = 2.23 --> % log_sigmay = 0.805
log_sigmax = log(sigmax);
log_sigmay = log(sigmay);


vrbl = [1:numdepvars];
% Now replicate the values in vrbl depending on the time value 
vrbl = repmat(vrbl,1,length(unique(time)));


%::::::::: GENERATE DATA ::::::::::::::::::::::::::::::::::::::::

bigtheta_true = [X0,log_sigmax,log_sigmay];   % store here all structural model parameters (log-scale)

model_param = {bigtheta_true,problem,time,numdepvars,vrbl};  


%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% create 'observations' by adding measurement noise a the noise-free
% trajectory
nobs = length(sampletime);
xhat_true = feval([problem, '_statemodel'],bigtheta_true,nobs,1); %noise free data (X)
yobs =  feval([problem, '_errormodel'],bigtheta_true,xhat_true,nobs,1); % noisy (observed) data Y

% plot(sampletime,yobs,'o-')
% xlabel('time')
% ylabel('Y')
% return
%::::::::: END OF DATA GENERATION ::::::::::::::::::::::::::::::::::::

                %X0   log_sigmax  log_sigmay
bigtheta_start = [0,     1.39,       1.39];  % mean value of the (log) starting parameter values
parmask        = [0 ,       1,         1 ];  % 1 for free parameters, 0 for constant/fixed ones
parbase = bigtheta_start;

% sample random starting values
theta_start_matrix = mvnrnd(bigtheta_start(2:3),[2 2],numattempts);  
bigtheta_start_matrix = [zeros(numattempts,1),theta_start_matrix];
all_theta_estimated = zeros(numattempts,sum(parmask));


% use the random starting parameter values
% estimate always from the same dataset
for ii=1:numattempts
    
    model_param{1} = bigtheta_start_matrix(ii,:);
    parbase = bigtheta_start_matrix(ii,:);
    % this is the SAEM-SL algorithm
    [THETAsaem,simsummaries] = saem_synlik(model_param,parmask,parbase,yobs,saem_numit,warmup,numsim);
    % save the parameter values obtained at the last iteration of SAEM-SL
    all_theta_estimated(ii,:) = THETAsaem(end,:);
    
    if ii==1
        figure
    end
    subplot(1,2,1)
    plot(exp(THETAsaem(:,1)))
    hline(sigmax)
    xlabel('\sigma_x')
    hold on
    subplot(1,2,2)
    plot(exp(THETAsaem(:,2)))
    hline(sigmay)
    xlabel('\sigma_y')
    hold on
end
saveas(gcf,'saem_trajectories')

save('all_theta_estimated','all_theta_estimated')
hold off






