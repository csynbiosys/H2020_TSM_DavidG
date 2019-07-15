%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TITLE: Toggle Switch model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%============================
% RESULTS PATHS RELATED DATA
%============================
inputs.pathd.results_folder='M1_'; % Folder to keep results (in Results) for a given problem                       
inputs.pathd.short_name='M1_RIdent';                       % To identify figures and reports for a given problem 
inputs.pathd.runident='full-par-bound';                         % [] Identifier required in order not to overwrite previous results

epcc_exps=1;


% Load experimental data. This data are derived, using the script
% ExtractionAllExperimentalData/ExtractionStructure_AllData_final.m, starting from data used in the paper
% Lugagne et al.
% The data derives from 26 experiments, 6 calibration and 20 control
% experiments.
% load('AllDataLugagne_Final.mat');
% load('AllDataLugagne_Final_minus_calib1_3.mat');
 load('AllDataLugagne_Final_21exp.mat');

%============================
% MODEL RELATED DATA
%============================
[model] = ToggleSwitch_load_model_M1();

%==================================
% EXPERIMENTAL SCHEME RELATED DATA
%==================================
 exps.n_exp=21;                                 % Number of experiments   
 
 % Store model details 
 inputs.model=model;
 inputs.model.exe_type='standard';

% Preprocess input
 AMIGO_Prep(inputs);

% Compute the steady state considering the initial theta guess, u_IPTG and
% u_aTc

Y0 = zeros(exps.n_exp,model.n_st);

for iexp=1:exps.n_exp
    Data.init_u{iexp}=[Data.Initial_IPTG{1,iexp}; Data.Initial_aTc{1,iexp}]; 
    Y0(iexp,:) = M1_Compute_SteadyState_OverNight(epcc_exps,inputs,iexp,model.par,Data.exp_data{1,iexp}(:,1)',Data.init_u{iexp});
end


  for iexp=1:exps.n_exp
    exps.exp_type{iexp} = 'fixed'; 
    exps.n_obs{iexp} = 2; 
     exps.obs_names{iexp} = char('LacI_M1','TetR_M1');
    exps.obs{iexp} = char('LacI_M1 = L_AU','TetR_M1 = T_AU');
    exps.t_in{iexp} = Data.t_con{1,iexp}(1,1); 
    exps.t_f{iexp} = Data.t_con{1,iexp}(1,end); 
    exps.n_s{iexp} = Data.n_samples{1,iexp}(1,1);
    exps.t_s{iexp} = Data.t_samples{1,iexp}(1,:); 
    exps.u_interp{iexp} = 'step';
    exps.t_con{iexp} = Data.t_con{1,iexp}(1,:);
    exps.n_steps{iexp} = length(exps.t_con{iexp})-1;
    exps.u{iexp} = Data.input{1,iexp};
    exps.data_type = 'real';
    exps.noise_type = 'homo';

    exps.exp_data{iexp} = Data.exp_data{1,iexp}';
    exps.error_data{iexp} = Data.standard_dev{1,iexp}';

    exps.exp_y0{iexp} = Y0(iexp,:);
    exps.exp_y0{iexp}([ 4 3]) = Data.exp_data{1,iexp}(:,1)';  
    %exps.init_u{iexp}=Data.init_u{iexp};
  end
    

% Store experiment details 
inputs.exps=exps;

%==================================
% UNKNOWNS RELATED DATA
%==================================

global_theta_min=[0.1386,5e-06,1e-3,40,1,2,2,0.0165,1e-05,0.001,0.0222,3.07e-4,2,2,1.32,1,0.001,0.001]; % verify Theta_T is correct 
global_theta_max=[0.1386,10,1000,40,100,4,4,0.0165,30,1000,0.0222,3.07e-2,4,4,1320,1000,0.1,0.1]; % Maximum allowed values for the parameters


% Parameters to be considered during sensitivity analysis
param_including_vector=[false,true,true,false,true,true,true,false,true,true,false,true,true,true,true,true,true,true];


% GLOBAL UNKNOWNS (SAME VALUE FOR ALL EXPERIMENTS)

inputs.PEsol.id_global_theta=model.par_names(param_including_vector,:); %'all'; %|User selected 
inputs.PEsol.global_theta_min=global_theta_min(param_including_vector); %0.5.*inputs.model.par(param_including_vector); %[0.1386,5e-06,1e-3,7000,1,2,2,0.0165,1e-05,0.001,700,3.07e-4,2,2,1.32,1,0.001,0.001]; % verify Theta_T is correct 
inputs.PEsol.global_theta_max=global_theta_max(param_including_vector); %2.*inputs.model.par(param_including_vector); %[0.1386,10,1000,7000,100,4,4,0.0165,30,1000,700,3.07e-2,4,4,1320,1000,0.1,0.1]; % Maximum allowed values for the parameters
inputs.PEsol.global_theta_guess=inputs.model.par(param_including_vector);

%==================================
% NUMERICAL METHODS RELATED DATA
%==================================

% SIMULATION
inputs.ivpsol.ivpsolver='cvodes';                     % [] IVP solver: C:'cvodes'; MATLAB:'ode15s'(default)|'ode45'|'ode113'            
inputs.ivpsol.senssolver='cvodes';                    % [] Sensitivities solver: 'cvodes' (C)                                                      
inputs.ivpsol.rtol=1e-8;                            % [] IVP solver integration tolerances
inputs.ivpsol.atol=1e-12; 


%==================================
% COST FUNCTION RELATED DATA
%==================================       
inputs.PEsol.PEcost_type='lsq';                       % 'lsq' (weighted least squares default) | 'llk' (log likelihood) | 'user_PEcost' 
inputs.PEsol.lsq_type='Q_expmax';
inputs.nlpsol.eSS.log_var = (1:14); 

% OPTIMIZATION
inputs.nlpsol.nlpsolver='eSS';                        % [] NLP solver: 

%==================================
% RIdent or GRank DATA
%==================================
inputs.rid.conf_ntrials=500;                          % [] Number of tria�ls for the robust confidence computation (default: 500)

inputs.nlpsol.eSS.local.solver = 'fminsearch';
inputs.nlpsol.eSS.local.finish = 'fminsearch';
inputs.nlpsol.eSS.local.nl2sol.maxfeval = 500;
inputs.nlpsol.eSS.maxeval = 5000;
inputs.nlpsol.eSS.maxtime = 500;
inputs.nlpsol.eSS.local.nl2sol.objrtol =  inputs.ivpsol.rtol;
inputs.nlpsol.eSS.local.nl2sol.tolrfun = 1e-12;

%==================================
% DISPLAY OF RESULTS
%==================================

inputs.plotd.plotlevel='min';                       % [] Display of figures: 'full'|'medium'(default)|'min' |'noplot' 

% Preprocess script
AMIGO_Prep(inputs)

AMIGO_RIdent(inputs)