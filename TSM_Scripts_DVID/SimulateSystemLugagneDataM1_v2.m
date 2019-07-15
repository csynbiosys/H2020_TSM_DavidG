epccOutputResultFileNameBase = 'Test1-1';
epcc_exps = 1;
global_theta_guess = [0.05938212,6.27986,92.18167,0.9259555,0.4403347,3.295859,0.2839974,4.49357,2.679409,...
    0.01641903,0.7912365,0.7186982,0.4464296,0.01388927];

LugDat = load('D:\PhD\GitHub\H2020_Varun\AMIGO2_R2018b_Varun\Examples\TSM\model_selection\M2\AllDataLugagne_Final.mat');
inexp = 5;

resultFileName = [strcat(epccOutputResultFileNameBase),'.dat'];
rng shuffle;
rngToGetSeed = rng;

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
short_name     = strcat('M1_Lugagne',int2str(epcc_exps));

model=ToggleSwitch_load_model_M1_DVID_v2();

exps.n_exp=0;

clear inputs;
inputs.model = model;
inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = 'initial_setup';
AMIGO_Prep(inputs);



y0 = M1_Compute_SteadyState_OverNight_DVID_v2(epcc_exps,inputs,1,global_theta_guess,...
    [LugDat.Data.exp_data{inexp}(1,1), LugDat.Data.exp_data{inexp}(2,1)],...
    [LugDat.Data.Initial_IPTG{inexp} LugDat.Data.Initial_aTc{inexp}]+1e-7);

duration = round(LugDat.Data.t_samples{inexp}(1,end));             
                                   
inducer = LugDat.Data.input{inexp}+1e-7;


clear newExps;
newExps.n_exp = 1;                                         % Number of experiments 
newExps.n_obs{1}=2;                                        % Number of observables per experiment                         
newExps.obs_names{1}=char('T_AU','L_AU');                 % Name of the observables per experiment    
newExps.obs{1}=char('LacI_AU=L_AU','TetR_AU=T_AU');           % Observation function
newExps.exp_y0{1}=y0;                                      % Initial condition for the experiment    

newExps.t_f{1}=duration;                                   % Experiment duration
newExps.n_s{1}=LugDat.Data.n_samples{inexp}(1);                             % Number of sampling times

sampl = round(LugDat.Data.t_samples{inexp}(1,10)-LugDat.Data.t_samples{inexp}(1,9));

newExps.t_s{1}=0:sampl:duration ;                              % Times of samples

newExps.u_interp{1}='step';                                % Interpolating function for the input

newExps.n_steps{1}=length(LugDat.Data.input{inexp}(1,:));                  % Number of steps in the input
newExps.u{1}= inducer;                                     % IPTG and aTc values for the input
newExps.t_con{1}=round(LugDat.Data.t_con{inexp}(1,:));                     % Switching times


clear inputs;
inputs.model = model;
inputs.exps = newExps;

inputs.exps.data_type='pseudo';
inputs.exps.noise_type='hetero_proportional';
inputs.exps.std_dev{1}=[0.0 0.0];

% SIMULATION
inputs.ivpsol.ivpsolver='cvodes';
inputs.ivpsol.senssolver='cvodes';
inputs.ivpsol.rtol=1.0D-13;
inputs.ivpsol.atol=1.0D-13;

inputs.plotd.plotlevel='medium';

inputs.pathd.results_folder = results_folder;                        
inputs.pathd.short_name     = short_name;
inputs.pathd.runident       = strcat('sim-',int2str(epcc_exps));

% simDa = AMIGO_SData(inputs);

simMo = AMIGO_SModel(inputs);









