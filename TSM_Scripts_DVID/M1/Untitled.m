clear;

epccOutputResultFileNameBase = 'Test1-1';
epcc_exps = 1;
global_theta_guess = [0.05938212,6.27986,92.18167,0.9259555,0.4403347,3.295859,0.2839974,4.49357,2.679409,...
    0.01641903,0.7912365,0.7186982,0.4464296,0.01388927];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit_to_InduciblePromoter_Step script - runs PE on the step class of inputs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resultFileName = [strcat(epccOutputResultFileNameBase),'.dat'];
rng shuffle;
rngToGetSeed = rng;

% Write the header information of the .dat file in which the results of
% PE (estimates, relative confidence intervals, residuals, relative
% residuals and the time required for computation) will be stored. 
fid = fopen(resultFileName,'w');
fprintf(fid,'HEADER DATE %s\n',datestr(datetime()));
fprintf(fid,'HEADER RANDSEED %d\n',rngToGetSeed.Seed);
fclose(fid);

startTime = datenum(now);

clear model;
clear exps;
clear best_global_theta;
clear pe_results;
clear pe_inputs;

% Specify folder name and short_name
results_folder = strcat('M1',datestr(now,'yyyy-mm-dd-HHMMSS'));
short_name     = strcat('M1_CIC_Random',int2str(epcc_exps));

% Read the model into the model variable
model=ToggleSwitch_load_model_M1_DVID_v2();

% Start with no experiments


% Compute the steady state considering the initial theta guess and 0 IPTG
%     y0 =  M1_compute_steady_state_DVID(global_theta_guess,[50.61466259329, 158.608748834354],[0 100]); %[ToggleSwitch_M1_Compute_SteadyState(model.par,[1000,400],1,0) 0];
% Define boundaries for the parameters (taken from scientific literature)
%     global_theta_min=[0.1386,5e-06,1e-3,7000,1,2,2,0.0165,1e-05,0.001,700,3.07e-4,2,2,1.32,1,0.001,0.001]; % verify Theta_T is correct 
%     global_theta_max=[0.1386,10,1000,7000,100,4,4,0.0165,30,1000,700,3.07e-2,4,4,1320,1000,0.1,0.1]; % Maximum allowed values for the parameters

global_theta_guess = global_theta_guess';

% Specify the parameters to be calibrated. 
% The selection depends on the identifiability analysis preceding the 
% comparison: parameters that are not identifiable will fixed to the best
% estimate available for them.
% In our case, all parameters are identifiable and for those parameters
% which are excluded, we know the values, and they are fixed.
%     param_including_vector=[false,true,true,false,true,true,true,false,true,true,false,true,true,true,true,true,true,true];

% Compile the model
clear inputs;
inputs.model = model;
inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = 'initial_setup';


AMIGO_Prep(inputs);

y0 = M1_Compute_SteadyState_OverNight_DVID_v2(epcc_exps,inputs,1,global_theta_guess,[23, 1400],[1 0]+1e-7);

iexp = 1;

inputs.exps.n_exp=1;
inputs.exps.n_obs{iexp}=2;                             % number of observables
inputs.exps.obs_names{iexp}=char('T_GFP','L_RFP');       % name of observables
inputs.exps.obs{iexp}=char('LacI_AU=L_RFP','TetR_AU=T_GFP'); 
inputs.exps.exp_y0{iexp}=y0;%M1_Compute_SteadyState_OverNight_DVID_v2(epcc_exps,inputs,1,global_theta_guess,[23, 1400],[1 0]+1e-7);

% Noise in the data 
inputs.exps.data_type='pseudo';
inputs.exps.noise_type='hetero_proportional';
inputs.exps.std_dev{iexp}=[0.0 0.0];


inputs.exps.u_interp{1}='step';                                % [] Stimuli definition: u_interp: 'sustained' |'step'|'linear'(default)|
inputs.exps.n_pulses{1}=7;
inputs.exps.t_f{1} = 20*60;
% inputs.exps.t_con{1}= 60.*[0 1.5 2 3.5 4 5.5 6 7.5 8 9.5 10 11.5 12 13.5 14 20];    % stimulus swithching times, experiment 1
inputs.exps.t_con{1}=[0 60*7 60*15 inputs.exps.t_f{1}];

inputs.exps.t_s{1}   = (0:5:(20*60));                         % sampling time is 5 minutes (5X60=300 seconds)
inputs.exps.n_s{1}   = 20*60/5+1; 
% inputs.exps.u_min{1}=[0.6; 0.01];  inputs.exps.u_max{1}=[0.6; 50];    
inputs.exps.u{1}=[0 1 0; 100 0 100]+1e-7;
inputs.exps.n_steps{1}=length(inputs.exps.t_con{1})-1;

global_theta_min = [0.03*0.1, 10*0.1, 30*0.1, 1, 0, 0,...
    0.1*0.1, 1*0.1, 30*0.1, 0.01, 0, 0,...
    (4e-2)*0.1, (1e-1)*0.1]; 
global_theta_max = [0.03*10, 10*10, 10*10, 100, 5, 5,...
    0.1*10, 1*10, 30*10, 1, 5, 5,...
    (4e-2)*10, (1e-1)*10];
global_theta_guess = global_theta_guess;

best_global_theta = global_theta_guess;
inputs.PEsol.id_global_theta=model.par_names;
inputs.PEsol.global_theta_guess=(best_global_theta);
inputs.PEsol.global_theta_max=global_theta_max;  % Maximum allowed values for the paramters
inputs.PEsol.global_theta_min=global_theta_min;  

% SIMULATION
inputs.ivpsol.ivpsolver='cvodes';                     % [] IVP solver: C:'cvodes'; MATLAB:'ode15s'(default)|'ode45'|'ode113'            
inputs.ivpsol.senssolver='cvodes';                    % [] Sensitivities solver: 'cvodes' (C)                                                      
inputs.ivpsol.rtol=1.0D-13;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0D-13; 


%==================================
% COST FUNCTION RELATED DATA
%==================================       
inputs.PEsol.PEcost_type='llk';                       % 'lsq' (weighted least squares default) | 'llk' (log likelihood) | 'user_PEcost' 
% inputs.PEsol.lsq_type='Q_expmax';
inputs.nlpsol.eSS.log_var =(1:14); 

% OPTIMIZATION
inputs.nlpsol.nlpsolver='eSS';                        % [] NLP solver: 
 
%==================================
% RIdent or GRank DATA
%==================================
inputs.rid.conf_ntrials=500;                          % [] Number of tria�ls for the robust confidence computation (default: 500)

inputs.nlpsol.eSS.local.solver = 'fmincon';
inputs.nlpsol.eSS.local.finish = 'fmincon';
inputs.nlpsol.eSS.local.nl2sol.maxfeval = 500;
inputs.nlpsol.eSS.maxeval = 5000;
inputs.nlpsol.eSS.maxtime = 500;
inputs.nlpsol.eSS.local.nl2sol.objrtol =  inputs.ivpsol.rtol;
inputs.nlpsol.eSS.local.nl2sol.tolrfun = 1e-13;

%==================================
% GRank number of samples
%================================== 
 inputs.rank.gr_samples=10000;                         % [] Number of samples for global sensitivities and global rank within LHS (default: 10000)    
 

%==================================
% DISPLAY OF RESULTS
%==================================

% inputs.plotd.plotlevel='noplot';                       % [] Display of figures: 'full'|'medium'(default)|'min' |'noplot' 

% Preprocess script
AMIGO_Prep(inputs)


results=AMIGO_SData(inputs);














