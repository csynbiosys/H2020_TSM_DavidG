epccOutputResultFileNameBase = 'Test1-1';
epcc_exps = 1;
global_theta_guess = [0.05938212,6.27986,92.18167,0.9259555,0.4403347,3.295859,0.2839974,4.49357,2.679409,...
    0.01641903,0.7912365,0.7186982,1,1,0.4464296,0.01388927];
 
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
model=ToggleSwitch_load_model_M1_DVID();

% Start with no experiments
exps.n_exp=0;

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

y0 = M1_Compute_SteadyState_OverNight_DVID(epcc_exps,inputs,1,global_theta_guess,[23, 1400],[1 0]+1e-7);


% Fixed parts of the experiment
duration = 24*60;               % Duration in of the experiment (minutes)
Stepd = 180;                                         % Duration of each step (minutes). Note that this value equals the response time, quantified in 80 mins for MPLac,r
inducer = [rand(1,round(duration/Stepd)); randi([0 100],1,round(duration/Stepd))]+1e-7;  % Extract 2 random integer values from the Uniform distribution in (0,1) ng/mL for IPTG and (0,100) for aTc to be used as low and high values in each cycle
%     inducer = [0; 100]+1e-7;
% inducer = [0.8 0.2; 20 80]+1e-7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a new experiment to simulate with the step input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear newExps;
newExps.n_exp = 1;                                         % Number of experiments 
newExps.n_obs{1}=2;                                        % Number of observables per experiment                         
newExps.obs_names{1}=char('T_AU','L_AU');                 % Name of the observables per experiment    
newExps.obs{1}=char('LacI_AU=L_AU','TetR_AU=T_AU');           % Observation function
newExps.exp_y0{1}=y0;                                      % Initial condition for the experiment    

newExps.t_f{1}=duration;                                   % Experiment duration
newExps.n_s{1}=duration/5 + 1;                             % Number of sampling times
newExps.t_s{1}=0:5:duration ;                              % Times of samples

newExps.u_interp{1}='step';                                % Interpolating function for the input
% newExps.u_interp{2}='step';                                % Interpolating function for the input
newExps.n_steps{1}=round(duration/Stepd);                  % Number of steps in the input
newExps.u{1}= inducer;                                     % IPTG and aTc values for the input
newExps.t_con{1}=[0:180:1440];                     % Switching times
%     newExps.t_con{1}=[0 3000];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mock the experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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



