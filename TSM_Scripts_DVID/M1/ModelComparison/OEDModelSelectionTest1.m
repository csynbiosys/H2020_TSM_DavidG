clear;


epccOutputResultFileNameBase = 'Test1-OEDms';
epcc_exps = 1;

%===========================
%RESULTS  PATHS RELATED DATA
%===========================

inputs.pathd.results_folder='ms1';  % Folder to keep results (in Results\)
inputs.pathd.short_name='ms1';      % To identify figures and reports


%======================
% MODEL RELATED DATA
%======================

model = ToggleSwitch_load_model_M1vs2_DVID();

inputs.model = model;


%==========================================
% Dynamic optimization problem formulation
%==========================================
global_theta_guess = [0.05938212,6.27986,92.18167,0.9259555,0.4403347,3.295859,...
    0.2839974,4.49357,2.679409,0.01641903,0.7912365,0.7186982,...
    0.4464296,...
    0.01388927,...
    6.01684739e-02, 6.77117550e+00, 1.86956182e+01, 9.64252299e-01,7.72604972e-01, 2.66915511e+00,...
    1.23423001e-02, 5.13166491e+00, 4.69759475e+00, 6.05248321e-02, 1.02453826e+00, 5.21761158e-01,...
    1.70391261e-02, 1.83532118e-01, ...
    4.79131949e-02, 1.24263999e-02];

inputs.model.par=global_theta_guess;

AMIGO_Prep(inputs);

y0 = M1vs2_Compute_SteadyState_OverNight_DVID(epcc_exps,inputs,1,global_theta_guess,...
    [28.510, 1363.193],...
    [1 0]+1e-7);

inputs.DOsol.y0=y0;                               %Initial conditions
inputs.DOsol.tf_type='fixed';                          %Process duration type: fixed or free
inputs.DOsol.tf_guess=24*60;                               %Process duration
% inputs.DOsol.t_f = 0:5:(24*60);

%COST FUNCTION
inputs.DOsol.DOcost_type='max';                        %Type of problem: max/min
inputs.DOsol.DOcost='abs(L_RFP-L_RFP2)';                              %Cost functional


%CVP (Control Vector Parameterization) DETAILS
inputs.DOsol.u_interp='stepf';                         %Control definition
                                                       %'sustained' |'stepf'|'step'|'linear'|
inputs.DOsol.n_steps=2;
inputs.DOsol.u_guess=[0.5 0.5; 50 50];% Initial guess for the input
inputs.DOsol.u_min= [0*ones(1,inputs.DOsol.n_steps); 0*ones(1,inputs.DOsol.n_steps)]+1e-7;
inputs.DOsol.u_max=[1*ones(1,inputs.DOsol.n_steps); 100*ones(1,inputs.DOsol.n_steps)];% Minimum and maximum value for the input
% inputs.DOsol.u_min=-0.3.*ones(1,inputs.DOsol.n_steps);
% inputs.DOsol.u_max=1.*ones(1,inputs.DOsol.n_steps);    % Minimum and maximum value for the input
inputs.DOsol.t_con=0:5/inputs.DOsol.n_steps:24*60;         % Input swithching times, including intial and
                                                       % final times
                                                       
%==================================
% NUMERICAL METHDOS RELATED DATA
%==================================

% SIMULATION
inputs.ivpsol.ivpsolver='cvodes';
inputs.ivpsol.senssolver='cvodes';

inputs.ivpsol.rtol=1.0D-13;
inputs.ivpsol.atol=1.0D-13;

%OPTIMIZATION
% inputs.nlpsol.nlpsolver='local_fmincon'; 
inputs.nlpsol.reopt='off'; 
inputs.nlpsol.nlpsolver='eSS';                    % [] NLP solver:

inputs.nlpsol.eSS.maxeval = 10000;               % Maximum number of cost function evaluations
inputs.nlpsol.eSS.maxtime = 30000;                  % Maximum computational time in seconds
inputs.nlpsol.eSS.local.solver = 'fmincon';       % Local solver- SQP
inputs.nlpsol.eSS.local.finish = 'fmincon';    % Local solver- Direct method

%================================
% CALL AMIGO2 from COMMAND LINE
%================================
% It is recommended to keep all inputs in a 'problem_file'.m.
% AMIGO2 DO task can be called as follows:

AMIGO_Prep(inputs);


AMIGO_DO(inputs);







