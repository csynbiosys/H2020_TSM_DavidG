%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TITLE: Toggle Switch model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epccOutputResultFileNameBase = 'Test1-1';
epcc_exps = 1;
global_theta_guess = [0.05938212,6.27986,92.18167,0.9259555,0.4403347,3.295859,0.2839974,4.49357,2.679409,...
    0.01641903,0.7912365,0.7186982,0.4464296,0.01388927];

%============================
% RESULTS PATHS RELATED DATA
%============================
inputs.pathd.results_folder='DVID_ToggleSwitch'; % Folder to keep results (in Results) for a given problem                       
inputs.pathd.short_name='DVIDTSM_ContourP';                       % To identify figures and reports for a given problem 
inputs.pathd.runident='run1';                         % [] Identifier required in order not to overwrite previous results
                                                     
%============================
% MODEL RELATED DATA
%============================
[model] = ToggleSwitch_load_model_M1_DVID_v2();

%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
 inputs.exps.n_exp=6;                                 % Number of experiments   
 
 % Store model details 
 inputs.model=model;
 inputs.model.exe_type='standard';
 AMIGO_Prep(inputs)
 
  for iexp= 1:inputs.exps.n_exp
     inputs.exps.n_obs{iexp}=2;                             % number of observables
     inputs.exps.obs_names{iexp}=char('T_GFP','L_RFP');       % name of observables
     inputs.exps.obs{iexp}=char('LacI_AU=L_RFP','TetR_AU=T_GFP'); 
     inputs.exps.exp_y0{iexp}=M1_Compute_SteadyState_OverNight_DVID_v2(epcc_exps,inputs,1,global_theta_guess,[23, 1400],[1 0]+1e-7);
    
     % Noise in the data 
     inputs.exps.data_type='pseudo';
     inputs.exps.noise_type='hetero_proportional';
     inputs.exps.std_dev{iexp}=[0.05 0.05];
     
 end
 
    
 for i3=1:2
    inputs.exps.u_interp{i3}='pulse-up';                                % [] Stimuli definition: u_interp: 'sustained' |'step'|'linear'(default)|
    inputs.exps.n_pulses{i3}=7;
    inputs.exps.t_f{i3} = 20*60;
    inputs.exps.t_con{i3}= 60.*[0 1.5 2 3.5 4 5.5 6 7.5 8 9.5 10 11.5 12 13.5 14 20];    % stimulus swithching times, experiment 1  
    inputs.exps.t_s{i3}   = (0:5:20*60);                         % sampling time is 5 minutes (5X60=300 seconds)
    inputs.exps.n_s{i3}   = length((0:5:20*60));                                % Number of sampling times

 end
    inputs.exps.u_min{1}=[0.01; 0.001]; inputs.exps.u_max{1}=[0.01; 100];    
    inputs.exps.u_min{2}=[0.6; 0.01];  inputs.exps.u_max{2}=[0.6; 50];    

    
for i2=3:6
    inputs.exps.u_interp{i2}='step';                                % [] Stimuli definition: u_interp: 'sustained' |'step'|'linear'(default)|
    inputs.exps.t_f{i2} = 20*60;
    inputs.exps.t_s{i2}   = (0:5:20*60);                         % sampling time is 5 minutes (5X60=300 seconds)
    inputs.exps.n_s{i2}   = length((0:5:20*60));                                % Number of sampling times
end
 

 inputs.exps.u{3}=[0 1 0; 100 0 100]+1e-7; inputs.exps.t_con{1,3}=[0 60*7 60*15 inputs.exps.t_f{1,3}];inputs.exps.n_steps{3}=length(inputs.exps.t_con{1,3})-1;
 inputs.exps.u{4}=[0 1;100 0]+1e-7; inputs.exps.t_con{1,4}=[0 60*7 inputs.exps.t_f{1,4}];inputs.exps.n_steps{4}=length(inputs.exps.t_con{1,4})-1;
 inputs.exps.u{5}=[1 0;0 100]+1e-7; inputs.exps.t_con{1,5}=[0 60*7 inputs.exps.t_f{1,5}];inputs.exps.n_steps{5}=length(inputs.exps.t_con{1,5})-1;
 inputs.exps.u{6}=[[1,1,0,1,1,0,1,1,0]+1e-7;[ 9,89,3,84,51,11,49,60,9]]; inputs.exps.t_con{1,6}=[0 60*1 60*2 60*5 60*7 60*12 60*15 60*18 inputs.exps.t_f{1,6}];inputs.exps.n_steps{6}=length(inputs.exps.t_con{1,6})-1;



%==================================
% UNKNOWNS RELATED DATA
%==================================

% GLOBAL UNKNOWNS (SAME VALUE FOR ALL EXPERIMENTS)

% global_theta_min=[0.1386,5e-06,1e-3,7000,1,2,2,0.0165,1e-05,0.001,700,3.07e-4,2,2,1.32,1,0.001,0.001]; % verify Theta_T is correct 
% global_theta_max=[0.1386,10,1000,7000,100,4,4,0.0165,30,1000,700,3.07e-2,4,4,1320,1000,0.1,0.1]; % Maximum allowed values for the parameters

% Excluding parameters which are fixed
% param_including_vector=[false,true,true,false,true,true,true,false,true,true,false,true,true,true,true,true,true,true];

% settings for PE
global_theta_min = [0.03*0.1, 10*0.1, 30*0.1, 0, 0, 0,...
    0.1*0.1, 1*0.1, 10*0.1, 0.01, 0, 0,...
    (5e-2)*0.1, (1e-1)*0.1]; 
global_theta_max = [0.03*10, 10*10, 10*10, 100, 5, 5,...
    0.1*10, 1*10, 30*10, 1, 5, 5,...
    (5e-2)*10, (1e-1)*10];
global_theta_guess = global_theta_guess;

best_global_theta = global_theta_guess;
inputs.PEsol.id_global_theta=model.par_names;
inputs.PEsol.global_theta_guess=(best_global_theta);
inputs.PEsol.global_theta_max=global_theta_max;  % Maximum allowed values for the paramters
inputs.PEsol.global_theta_min=global_theta_min;  

%==================================
% NUMERICAL METHODS RELATED DATA

%==================================

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

inputs.plotd.plotlevel='noplot';                       % [] Display of figures: 'full'|'medium'(default)|'min' |'noplot' 

% Preprocess script
AMIGO_Prep(inputs)

% generate pseudo data
results=AMIGO_SData(inputs);

for iexp=1:inputs.exps.n_exp
    inputs.exps.exp_data{1,iexp}=results.sim.exp_data{1,iexp};
    inputs.exps.error_data{1,iexp}=results.sim.error_data{1,iexp};  
end

inputs.plotd.nx_contour=200;                          % Number of points for plotting the contours x and y direction
inputs.plotd.ny_contour=200;                          % ADVISE: >50

AMIGO_ContourP(inputs)