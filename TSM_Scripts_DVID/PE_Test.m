epccOutputResultFileNameBase = 'Test1-1';
epcc_exps = 1;
global_theta_guess = [0.05938212,6.27986,92.18167,0.9259555,0.4403347,3.295859,0.2839974,4.49357,2.679409,...
    0.01641903,0.7912365,0.7186982,0.4464296,0.01388927];

% Specify folder name and short_name
results_folder = strcat('M1Fit',datestr(now,'yyyy-mm-dd-HHMMSS'));
short_name     = strcat('M1F',int2str(epcc_exps));

% Load data
LugDat = load('D:\PhD\GitHub\H2020_Varun\AMIGO2_R2018b\Examples\TSM\model_selection\M2\AllDataLugagne_Final.mat');
inexp = 14;

% Load model
model=ToggleSwitch_load_model_M1_DVID_v2();

% Initial guesses for theta 
global_theta_min = [0.03*0.1, 10*0.1, 30*0.1, 1, 0, 0,...
    0.1*0.1, 1*0.1, 30*0.1, 0.01, 0, 0,...
    (4e-2)*0.1, (1e-1)*0.1]; 
global_theta_max = [0.03*10, 10*10, 10*10, 100, 5, 5,...
    0.1*10, 1*10, 30*10, 1, 5, 5,...
    (4e-2)*10, (1e-1)*10];
global_theta_guess = global_theta_guess;


% Compile the model
clear inputs;
inputs.model = model;
inputs.model.par = global_theta_guess;

inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = 'initial_setup';
AMIGO_Prep(inputs);

% Compute steady state
y0 = M1_Compute_SteadyState_OverNight_DVID_v2(epcc_exps,inputs,1,global_theta_guess,...
    [LugDat.Data.exp_data{inexp}(1,1), LugDat.Data.exp_data{inexp}(2,1)],...
    [LugDat.Data.Initial_IPTG{inexp} LugDat.Data.Initial_aTc{inexp}]+1e-7);

clear exps;

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


best_global_theta = global_theta_guess; 

% Compile the model
clear inputs;
inputs.model = model;
inputs.model.par = best_global_theta;
inputs.exps  = exps;

inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = strcat('pe-',int2str(epcc_exps));

AMIGO_Prep(inputs);


% GLOBAL UNKNOWNS (SAME VALUE FOR ALL EXPERMENTS)
inputs.PEsol.id_global_theta=model.par_names;
inputs.PEsol.global_theta_guess=(best_global_theta);
inputs.PEsol.global_theta_max=global_theta_max;  % Maximum allowed values for the paramters
inputs.PEsol.global_theta_min=global_theta_min;  % Minimum allowed values for the parameters


 %COST FUNCTION RELATED DATA
inputs.PEsol.PEcost_type='llk';                       % 'lsq' (weighted least squares default) | 'llk' (log likelihood) | 'user_PEcost'
inputs.PEsol.llk_type='homo_var';    

%SIMULATION
inputs.ivpsol.ivpsolver='cvodes';
inputs.ivpsol.senssolver='cvodes';
inputs.ivpsol.rtol=1.0D-13;
inputs.ivpsol.atol=1.0D-13;

% OPTIMIZATION
% inputs.nlpsol.nlpsolver='hyb_de_nl2sol';    % In this case the problem will be solved with a
%                                             % sequential hybrid - DE + nl2sol
% inputs.nlpsol.DE.NP=100;                    % Size of population - recommended 10*N variables
% inputs.nlpsol.DE.itermax=200;               % Maximum number of iterations - 100*200 function
%                                             % evaluations
% inputs.nlpsol.DE.strategy=7;                % Strategy in DE 2:DE/rand/1/exp
% inputs.nlpsol.DE.cvarmax=1e-4;              % Maximum variance for the population


inputs.nlpsol.nlpsolver='eSS';
inputs.nlpsol.eSS.maxeval = 200000;
inputs.nlpsol.eSS.maxtime = 5000;
inputs.nlpsol.eSS.local.solver = 'fmincon'; 
inputs.nlpsol.eSS.local.finish = 'fmincon'; 
inputs.rid.conf_ntrials=500;


pe_start = now;
pe_inputs = inputs;
results = AMIGO_PE(inputs);
pe_results = results;
pe_end = now;

%Save the best theta
best_global_theta = results.fit.thetabest;

%Write results to the output file
fid = fopen(resultFileName,'a');
used_par_names = model.par_names;

for j=1:size(used_par_names,1)
    fprintf(fid,'PARAM_FIT %s %f\n', used_par_names(j,:), results.fit.thetabest(j));
end

%Time in seconds
fprintf(fid,'PE_TIME %.1f\n', (pe_end-pe_start)*24*60*60);
fclose(fid);

save(strcat(epccOutputResultFileNameBase,'PE2.mat'),'pe_results','exps','pe_inputs','best_global_theta');






























