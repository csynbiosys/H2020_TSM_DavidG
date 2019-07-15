
clear;

epccOutputResultFileNameBase = 'Test1-1';
epcc_exps = 1;
global_theta_guess = [0.05938212,6.27986,92.18167,0.9259555,0.4403347,3.295859,0.2839974,4.49357,2.679409,...
    0.01641903,0.7912365,0.7186982,0.4464296,0.01388927];

inputs.pathd.results_folder='M1LSens';  % Folder to keep results (in Results\)
inputs.pathd.short_name='M1LSens';      % Label to identify figures and reports

LugDat = load('D:\PhD\GitHub\H2020_Varun\AMIGO2_R2018b_Varun\Examples\TSM\model_selection\M2\AllDataLugagne_Final.mat');
inexp = 14;

model = ToggleSwitch_load_model_M1_DVID_v2();

inputs.model = model;

inputs.model.par=global_theta_guess;

AMIGO_Prep(inputs);


y0 = M1_Compute_SteadyState_OverNight_DVID_v2(epcc_exps,inputs,1,global_theta_guess,...
    [26 1400],...
    [1 0]+1e-7);

inputs.exps.n_exp=1;                          % Number of experiments
duration = 24*60;
Stepd = 720;  
inducer = [0 1; 100 0]+1e-7;
% EXPERIMENT 1

inputs.exps.exp_y0{1}=y0;        % Initial conditions
inputs.exps.t_f{1}=duration;                       % Experiments duration

inputs.exps.n_obs{1}=2;                       % Number of observables
inputs.exps.obs_names{1}=char('RFP','GFP'); % Names of the observables
inputs.exps.obs{1}=char('RFP = L_RFP','GFP = T_GFP');

inputs.exps.u_interp{1}='step';          % Stimuli definition for experiment 1
inputs.exps.t_con{1}=[0 120];                 % Input swithching times including:
inputs.exps.n_s{1}=duration/5 + 1;                             % Number of sampling times
inputs.exps.t_s{1}=0:5:duration ;
inputs.exps.n_steps{1}=round(duration/Stepd);                  % Number of steps in the input
inputs.exps.u{1}= inducer;                                     % IPTG and aTc values for the input
inputs.exps.t_con{1}=[0:720:1440];                     % Switching times

%============================================
% PARAMETERS TO BE CONSIDERED IN THE ANALYSIS
%============================================

global_theta_min = [0.03*0.1, 10*0.1, 30*0.1, 0, 0, 0,...
    0.1*0.1, 1*0.1, 10*0.1, 0.01, 0, 0,...
    (5e-2)*0.1, (1e-1)*0.1]; 
global_theta_max = [0.03*10, 10*10, 10*10, 100, 5, 5,...
    0.1*10, 1*10, 30*10, 1, 5, 5,...
    (5e-2)*10, (1e-1)*10];
best_global_theta = global_theta_guess;

inputs.PEsol.id_global_theta=model.par_names;
inputs.PEsol.global_theta_guess=(best_global_theta);
inputs.PEsol.global_theta_max=global_theta_max;  % Maximum allowed values for the paramters
inputs.PEsol.global_theta_min=global_theta_min;  % Minimum allowed values for the parameters


inputs.rank.gr_samples=10000;

AMIGO_Prep(inputs);

AMIGO_GRank(inputs);







%%%%%%%%%%%%%%%%%%% Real Experimental data
clear;

epccOutputResultFileNameBase = 'Test1-1';
epcc_exps = 1;
global_theta_guess = [0.05938212,6.27986,92.18167,0.9259555,0.4403347,3.295859,0.2839974,4.49357,2.679409,...
    0.01641903,0.7912365,0.7186982,0.4464296,0.01388927];

inputs.pathd.results_folder='M1LSens';  % Folder to keep results (in Results\)
inputs.pathd.short_name='M1LSens';      % Label to identify figures and reports

LugDat = load('D:\PhD\GitHub\H2020_Varun\AMIGO2_R2018b_Varun\Examples\TSM\model_selection\M2\AllDataLugagne_Final.mat');
inexp = 14;

model = ToggleSwitch_load_model_M1_DVID_v2();

inputs.model = model;

inputs.model.par=global_theta_guess;

AMIGO_Prep(inputs);

%============================================
% EXPERIMENTAL SCHEME (SIMULATION CONDITIONS)
%============================================
y0 = M1_Compute_SteadyState_OverNight_DVID_v2(epcc_exps,inputs,1,global_theta_guess,...
    [LugDat.Data.exp_data{inexp}(1,1), LugDat.Data.exp_data{inexp}(2,1)],...
    [LugDat.Data.Initial_IPTG{inexp} LugDat.Data.Initial_aTc{inexp}]+1e-7);

inputs.exps.n_exp=1;                          % Number of experiments



duration = round(LugDat.Data.t_samples{inexp}(1,end));                                                
inducer = LugDat.Data.input{inexp}+1e-7;

% Define the exps structure, containing the experimental data to fit
exps.n_exp = length(inexp);
exp_indexData = inexp;
iexp = 1;
exps.exp_type{iexp} = 'fixed'; 

exps.n_obs{iexp} = 2; 
exps.obs_names{iexp} = char('RFP','GFP');
exps.obs{iexp} = char('RFP = L_RFP','GFP = T_GFP');

exps.t_f{iexp} = round(LugDat.Data.t_samples{inexp}(1,end));
exps.n_s{iexp} = LugDat.Data.n_samples{inexp}(1);
sampl = round(LugDat.Data.t_samples{inexp}(1,10)-LugDat.Data.t_samples{inexp}(1,9));
exps.t_s{iexp} = 0:sampl:duration; 

exps.u_interp{iexp} = 'step';
exps.t_con{iexp} = round(LugDat.Data.t_con{inexp}(1,:)); 
exps.n_steps{iexp} = length(LugDat.Data.input{inexp}(1,:));
exps.u{iexp} = inducer;
exps.data_type = 'real';
exps.noise_type = 'homo';

exps.exp_data{iexp} = LugDat.Data.exp_data{inexp}';
exps.error_data{iexp} = LugDat.Data.standard_dev{inexp}';
exps.exp_y0{iexp} = y0(1,:);

inputs.exps  = exps;

%============================================
% PARAMETERS TO BE CONSIDERED IN THE ANALYSIS
%============================================

global_theta_min = [0.03*0.1, 10*0.1, 30*0.1, 0, 0, 0,...
    0.1*0.1, 1*0.1, 10*0.1, 0.01, 0, 0,...
    (5e-2)*0.1, (1e-1)*0.1]; 
global_theta_max = [0.03*10, 10*10, 10*10, 100, 5, 5,...
    0.1*10, 1*10, 30*10, 1, 5, 5,...
    (5e-2)*10, (1e-1)*10];
best_global_theta = global_theta_guess;

inputs.PEsol.id_global_theta=model.par_names;
inputs.PEsol.global_theta_guess=(best_global_theta);
inputs.PEsol.global_theta_max=global_theta_max;  % Maximum allowed values for the paramters
inputs.PEsol.global_theta_min=global_theta_min;  % Minimum allowed values for the parameters


inputs.rank.gr_samples=5000;

AMIGO_Prep(inputs);

AMIGO_GRank(inputs);



