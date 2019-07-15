
clear;

epccOutputResultFileNameBase = 'Test1-OED';
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
    [LugDat.Data.exp_data{inexp}(1,1), LugDat.Data.exp_data{inexp}(2,1)],...
    [LugDat.Data.Initial_IPTG{inexp} LugDat.Data.Initial_aTc{inexp}]+1e-7);


inputs.exps.n_exp=2;                             % Total number of experiments
                                                  %(available + experiments to be designed)
inputs.exps.exp_type{1}='fixed';                 % Indicates if the the experiment should
inputs.exps.exp_type{2}='od';

for iexp=1:inputs.exps.n_exp
    % OBSEVABLES DEFINITION
    inputs.exps.n_obs{iexp}=2;                       % Number of observed quantities per experiment
    inputs.exps.obs_names{iexp}=char('RFP','GFP'); % Name of the observed quantities per experiment
    inputs.exps.obs{iexp}=char('RFP = L_RFP','GFP = T_GFP');
    inputs.exps.exp_y0{iexp}=zeros(1,inputs.model.n_st);% Initial conditions for each experiment
 end

duration = round(LugDat.Data.t_samples{inexp}(1,end));                                                
inducer = LugDat.Data.input{inexp}+1e-7;

% Define the exps structure, containing the experimental data to fit
exp_indexData = inexp;

inputs.exps.t_f{1} = round(LugDat.Data.t_samples{inexp}(1,end));
inputs.exps.n_s{1} = LugDat.Data.n_samples{inexp}(1);
sampl = round(LugDat.Data.t_samples{inexp}(1,10)-LugDat.Data.t_samples{inexp}(1,9));
inputs.exps.t_s{1} = 0:sampl:duration; 

inputs.exps.u_interp{1} = 'step';
inputs.exps.t_con{1} = round(LugDat.Data.t_con{inexp}(1,:)); 
inputs.exps.n_steps{1} = length(LugDat.Data.input{inexp}(1,:));
inputs.exps.u{1} = inducer;
inputs.exps.data_type = 'real';
inputs.exps.noise_type = 'homo';

inputs.exps.exp_data{1} = LugDat.Data.exp_data{inexp}';
inputs.exps.error_data{1} = LugDat.Data.standard_dev{inexp}';
inputs.exps.exp_y0{1} = y0(1,:);



%
%  INPUTS FOR THE EXPERIMENT TO BE OPTIMALLY DESIGNED
%

 inputs.exps.u_type{2}='od';                       % Stimulation: 'fixed' | 'od' (to be designed)
 inputs.exps.u_interp{2}='step';                   % Stimuli definition for experiment 3:
                                                   % OPTIONS:u_interp: 'sustained' |'step'|
                                                   %        'linear'(default)|'pulse-up'|'pulse-down'
 inputs.exps.n_steps{2}=8;                         % Number of pulses _|-|_|-|_
 inputs.exps.u_min{2}= [0*ones(1,inputs.exps.n_steps{2}); 0*ones(1,inputs.exps.n_steps{2})]+1e-7;
 inputs.exps.u_max{2}=[1*ones(1,inputs.exps.n_steps{2}); 100*ones(1,inputs.exps.n_steps{2})];% Minimum and maximum value for the input
 inputs.exps.tf_type{2}='fixed';                   % [] Experiment duration:'fixed'(def.)|'od'(to be designed)
 inputs.exps.t_f{2}=24*60;                           % Experiment duration
 inputs.exps.ts_type{2}='fixed';                   % [] Sampling times:'fixed'(def.)| 'od'(to be designed)
 inputs.exps.n_s{2}=duration/5 + 1;
 inputs.exps.std_dev{2}=0.05;                       % Standard deviation of the noise for each
 inputs.exps.t_s{1}=0:5:duration ;                                                  % experiment: Ex: 0.05 <=> 5%

 
%======================================
% PARAMETERS TO BE CONSIDERED FOR OED
%======================================


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


 %COST FUNCTION RELATED DATA
% inputs.PEsol.PEcost_type='llk';                       % 'lsq' (weighted least squares default) | 'llk' (log likelihood) | 'user_PEcost'
% inputs.PEsol.llk_type='homo_var';    

%SIMULATION
% inputs.ivpsol.ivpsolver='cvodes';
% inputs.ivpsol.senssolver='cvodes';
% inputs.ivpsol.rtol=1.0D-13;
% inputs.ivpsol.atol=1.0D-13;

%==================================
% COST FUNCTION RELATED DATA
%==================================

inputs.exps.noise_type='homo_var';                % Experimental noise: 'homo' |'homo_var'| 'hetero'
inputs.OEDsol.OEDcost_type='Dopt';              % FIM criterium: 'Dopt'|'Eopt'|'Aopt'|'Emod'|'DoverE'


%==================================
% NUMERICAL METHODS RELATED DATA
%==================================

% OPTIMIZATION

inputs.nlpsol.nlpsolver='eSS';                    % [] NLP solver:

inputs.nlpsol.eSS.maxeval = 100;               % Maximum number of cost function evaluations
inputs.nlpsol.eSS.maxtime = 300;                  % Maximum computational time in seconds
inputs.nlpsol.eSS.local.solver = 'fmincon';       % Local solver- SQP
inputs.nlpsol.eSS.local.finish = 'fmincon';    % Local solver- Direct method


%====================================================
% CALL AMIGO2 from COMMAND LINE - STEP-WISE DESIGN
%====================================================
% It is recommended to keep all inputs in a 'problem_file'.m.
% AMIGO2 OED task can be called as follows:
% AMIGO_OED('problem_file','run_ident') or AMIGO_OED(inputs)


AMIGO_Prep(inputs);


AMIGO_OED(inputs);


























