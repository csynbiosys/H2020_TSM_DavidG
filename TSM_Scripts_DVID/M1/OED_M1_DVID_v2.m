function [out]=OED_M1_DVID(epccOutputResultFileNameBase,epccNumLoops,stepd,epcc_exps,global_theta_guess)
    
    resultFileName = [strcat(epccOutputResultFileNameBase),'.dat'];
    rng shuffle;
    rngToGetSeed = rng;

    % Write the header information of the .dat file in which the results of
    % PE (estimates, relative confidence intervals, residuals, relative
    % residuals and the time required for computation) will be stored.
    fid = fopen(resultFileName,'w');
    fprintf(fid,'HEADER DATE %s\n', datestr(datetime()));
    fprintf(fid,'HEADER RANDSEED %d\n',rngToGetSeed.Seed);
    fclose(fid);

    startTime = datenum(now);

    clear model;
    clear exps;
    clear best_global_theta_log;
    clear pe_results;
    clear ode_results;

    results_folder = strcat('M1_OED',datestr(now,'yyyy-mm-dd-HHMMSS'));
    short_name     = strcat('M1_OED',int2str(epcc_exps));

    model=ToggleSwitch_load_model_M1_DVID_v2();
    
    % We start with no experiments
    exps.n_exp=0;

    % Define boundaries for the parameters 
    global_theta_min=[0.03*0.1, 10*0.1, 10*0.1, 0, 0, 0,...
    0.1*0.1, 1*0.1, 30*0.1, 0.01, 0, 0,...
    (5e-2)*0.1, (1e-1)*0.1]; 

    global_theta_max=[0.03*10, 10*10, 10*10, 100, 5, 5,...
    0.1*10, 1*10, 30*10, 1, 5, 5,...
    (5e-2)*10, (1e-1)*10];


    global_theta_guess = global_theta_guess';
    
    % Compile the model
    clear inputs;
    inputs.model = model;
    inputs.pathd.results_folder = results_folder;
    inputs.pathd.short_name     = short_name;
    inputs.pathd.runident       = 'initial_setup';
    AMIGO_Prep(inputs);
    
    % Fixed parts of the experiment
    totalDuration = 24*60;               % Duration in of the experiment (minutes)
    numLoops = epccNumLoops;             % Number of OID loops
    duration = totalDuration/numLoops;   % Duration of each loop (in our case the number is 1)
    stepDuration = stepd;
    
    
    for i=1:numLoops
    
    %Compute the steady state considering the initial theta guess and 0 IPTG
    
    y0 = M1_Compute_SteadyState_OverNight_DVID_v2(epcc_exps,inputs,1,global_theta_guess,[23, 1400],[1 0]+1e-7);
    
    % --------------------- This is for On-Line, will take care of it
    % later
    
    % Need to determine the starting state of the next part of the
    % experiment we wish to design the input for. If there are multiple loops,
    % We get this state by simulating the model using the current best theta
    % for the duration for which we have designed input.
%     if exps.n_exp == 0
%         oid_y0 = [y0]; 
%         best_global_theta = global_theta_guess;
%     else
%         % Simulate the experiment without noise to find end state
%         clear inputs;
%         inputs.model = model;
%         inputs.model.par = best_global_theta;
%         inputs.exps = exps;
%         
%         inputs.plotd.plotlevel='noplot';
%         inputs.pathd.results_folder = results_folder;
%         inputs.pathd.short_name     = short_name;
%         inputs.pathd.runident       = strcat('sim-',int2str(i));
%         
%         sim = AMIGO_SData(inputs);
%         
%         oid_y0 = [sim.sim.states{1}(end,:)];
%         
%     end

    oid_y0 = y0; 
    best_global_theta = global_theta_guess;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimal experiment design
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear inputs;
    inputs.model = model;
    inputs.exps  = exps;
%     format long g
    
    
    % Add new experiment that is to be designed
    inputs.exps.n_exp = inputs.exps.n_exp + 1;                      % Number of experiments
    iexp = inputs.exps.n_exp;                                       % Index of the experiment
    inputs.exps.exp_type{iexp}='od';                                % Specify the type of experiment: 'od' optimally designed
    inputs.exps.n_obs{iexp}=2;                                      % Number of observables in the experiment
    inputs.exps.obs_names{iexp}=char('RFP','GFP');               % Name of the observables in the experiment
    inputs.exps.obs{iexp}=char('RFP = L_RFP','GFP = T_GFP');          % Observation function
    
    
    % Fixed parts of the experiment
    inputs.exps.exp_y0{iexp}=oid_y0;                                % Initial conditions
    inputs.exps.t_f{iexp}=duration;                                 % Duration of the experiment (minutes)
    inputs.exps.n_s{iexp}=duration/5+1;                             % Number of sampling times - sample every 5 min
    
    % OED of the input
    inputs.exps.u_type{iexp}='od';
    inputs.exps.u_interp{iexp}='stepf';                             % Stimuli definition for experiment: 'stepf' steps of constant duration
    inputs.exps.n_steps{iexp}=round(duration/stepDuration);         % Number of steps in the input
    inputs.exps.t_con{iexp}=0:stepDuration:(duration);            % Switching times
    
    inputs.exps.u_min{iexp}= [0*ones(1,inputs.exps.n_steps{iexp}); 0*ones(1,inputs.exps.n_steps{iexp})]+1e-7;
    inputs.exps.u_max{iexp}=[1*ones(1,inputs.exps.n_steps{iexp}); 100*ones(1,inputs.exps.n_steps{iexp})];
    
    
    
    inputs.PEsol.id_global_theta=model.par_names;
    inputs.PEsol.global_theta_guess=transpose(global_theta_guess);
    inputs.PEsol.global_theta_max=global_theta_max;  % Maximum allowed values for the parameters
    inputs.PEsol.global_theta_min=global_theta_min;  % Minimum allowed values for the parameters
    
    
    inputs.exps.noise_type='hetero_proportional';           % Experimental noise type: Homoscedastic: 'homo'|'homo_var'(default)
    inputs.exps.std_dev{iexp}=[0.05 0.05];
    inputs.OEDsol.OEDcost_type='Dopt';
    
    
    % SIMULATION
%     inputs.ivpsol.ivpsolver='cvodes';                     % [] IVP solver: 'cvodes'(default, C)|'ode15s' (default, MATLAB, sbml)|'ode113'|'ode45'
%     inputs.ivpsol.senssolver='cvodes';                    % [] Sensitivities solver: 'cvodes'(default, C)| 'sensmat'(matlab)|'fdsens2'|'fdsens5'
%     inputs.ivpsol.rtol=1.0D-13;                            % [] IVP solver integration tolerances
%     inputs.ivpsol.atol=1.0D-13;
    
    % OPTIMIZATION
    %oidDuration=600;
    inputs.nlpsol.nlpsolver='eSS';
    inputs.nlpsol.eSS.maxeval = 200000;
    inputs.nlpsol.eSS.maxtime = 5000;
    inputs.nlpsol.eSS.local.solver = 'fmincon'; 
    inputs.nlpsol.eSS.local.finish = 'fmincon';
    
    inputs.nlpsol.eSS.local.nl2sol.maxiter  =     500;     % max number of iteration
    inputs.nlpsol.eSS.local.nl2sol.maxfeval =     500;     % max number of function evaluation
%     inputs.nlpsol.eSS.log_var=1:inputs.exps.n_steps{iexp};
%     inputs.plotd.plotlevel='min';
    
    inputs.pathd.results_folder = results_folder;
    inputs.pathd.short_name     = short_name;
    inputs.pathd.runident       = strcat('oed-',int2str(i));
        
    oed_start = now;
    
    results = AMIGO_OED(inputs);
    oed_results{i} = results;
    oed_end = now;
    
    results.plotd.plotlevel = 'noplot';
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create a new experiment to simulate with the merged input
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear newExps;
    newExps.n_exp = 1;
    newExps.n_obs{iexp}=2;                                      % Number of observables in the experiment
    newExps.obs_names{iexp}=char('RFP','GFP');               % Name of the observables in the experiment
    newExps.obs{iexp}=char('RFP = L_RFP','GFP = T_GFP');          % Observation function
    newExps.exp_y0{iexp}= [y0];
    
    newExps.t_f{iexp}=i*duration;
    newExps.n_s{iexp}=(i*duration)/5 + 1;
    newExps.t_s{iexp}=0:5:(i*duration);
    
    newExps.u_interp{iexp}='step';
    newExps.n_steps{iexp}=(i*duration)/stepDuration;
    newExps.t_con{iexp}=0:stepDuration:(i*duration);
    
    % Merge the input signal
    
    if exps.n_exp == 0
        newExps.u{iexp}=results.oed.u{results.oed.n_exp};
    else
        newExps.u{iexp}=[exps.u{iexp} results.oed.u{results.oed.n_exp}];
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mock an experiment
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear inputs;
    inputs.model = model;
    inputs.exps = newExps;
    
    inputs.exps.data_type='pseudo';
    inputs.exps.noise_type='hetero_proportional';
    inputs.exps.std_dev{1}=[0.05 0.05];
    
    inputs.plotd.plotlevel='noplot';
    
    inputs.pathd.results_folder = results_folder;
    inputs.pathd.short_name     = short_name;
    inputs.pathd.runident       = strcat('sim-',int2str(epcc_exps),'-',int2str(i));
    
    sim = AMIGO_SData(inputs);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now we need to add this experiment output to newExps and copy
    % to exps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    newExps.exp_type{1}='fixed';
    if exps.n_exp == 0
        newExps.exp_data{1}   = sim.sim.exp_data{1};
        newExps.error_data{1} = sim.sim.error_data{1};
    else
        newExps.exp_data{1}   = [ exps.exp_data{1}   ; sim.sim.exp_data{1}((size(exps.exp_data{1},1)+1):end)];
        newExps.error_data{1} = [ exps.error_data{1} ; sim.sim.error_data{1}((size(exps.exp_data{1},1)+1):end)];
    end
    
    newExps.data_type='real';                                     % Type of data: 'pseudo'|'pseudo_pos'|'real'
    newExps.noise_type='hetero_proportional';                                % Experimental noise type: Homoscedastic: 'homo'|'homo_var'(default)
    
    exps = newExps;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameter estimation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear inputs;
    inputs.model = model;
    inputs.exps  = exps;
    
    inputs.pathd.results_folder = results_folder;
    inputs.pathd.short_name     = short_name;
    inputs.pathd.runident       = strcat('pe-',int2str(epcc_exps),'-',int2str(i));
    
    % GLOBAL UNKNOWNS (SAME VALUE FOR ALL EXPERIMENTS)
    inputs.PEsol.id_global_theta=model.par_names;
    inputs.PEsol.global_theta_guess=transpose(global_theta_guess);
    inputs.PEsol.global_theta_max=global_theta_max;  % Maximum allowed values for the paramters
    inputs.PEsol.global_theta_min=global_theta_min;  % Minimum allowed values for the parameters
    
    % COST FUNCTION RELATED DATA
    inputs.PEsol.PEcost_type='llk';                       % 'lsq' (weighted least squares default) | 'llk' (log likelihood) | 'user_PEcost'
    inputs.PEsol.llk_type='homo_var'; 
    
    % SIMULATION
    inputs.ivpsol.ivpsolver='cvodes';
    inputs.ivpsol.senssolver='cvodes';
    inputs.ivpsol.rtol=1.0D-13;
    inputs.ivpsol.atol=1.0D-13;
    
    % OPTIMIZATION
    inputs.nlpsol.nlpsolver='eSS';
    inputs.nlpsol.eSS.maxeval = 2000;%200000;
    inputs.nlpsol.eSS.maxtime = 2000;
    inputs.nlpsol.eSS.local.solver = 'fmincon';  
    inputs.nlpsol.eSS.local.finish = 'fmincon';  
    inputs.rid.conf_ntrials=500;
    
    inputs.plotd.plotlevel='min';
    
    pe_start = now;
    results = AMIGO_PE(inputs);
    pe_inputs{i} = inputs;
    pe_results{i} = results;
    pe_end= now;
    
    results.plotd.plotlevel = 'min';
    
     %Save the best theta
    best_global_theta=results.fit.thetabest;
    
    %Write some results to the output file
    fid = fopen(resultFileName,'a');
    used_par_names = model.par_names;
    
    for j=1:size(used_par_names,1)
        fprintf(fid,'ITERATION %d PARAM_FIT %s %f\n', i, used_par_names(j,:), results.fit.thetabest(j));
        if isfield(results.fit,'rel_conf_interval')
            fprintf(fid,'ITERATION %d REL_CONF %s %f\n',  i, used_par_names(j,:), results.fit.rel_conf_interval(j));
        end
        if isfield(results.fit,'residuals')
            fprintf(fid,'ITERATION %d RESIDUAL %s %f\n', i, used_par_names(j,:), results.fit.residuals{1}(j));
        end
        if isfield(results.fit,'rel_residuals')
            fprintf(fid,'ITERATION %d REL_RESIDUAL %s %f\n', i, used_par_names(j,:), results.fit.rel_residuals{1}(j));
        end
    end
    %Time in seconds
    fprintf(fid,'ITERATION %d OED_TIME %.1f\n', i, (oed_end-oed_start)*24*60*60);
    fprintf(fid,'ITERATION %d PE_TIME %.1f\n',  i, (pe_end-pe_start)*24*60*60);
    fclose(fid);
    

    end
    
for i=1:8

    duration = i*3*60;  % Duration in minutes

    clear inputs;
    inputs.model = model;
    inputs.exps  = exps;

    % Reduce the input to a smaller set of values
    inputs.exps.t_f{1}          = duration;                % Experiment duration
    inputs.exps.n_s{1}          = duration/5 + 1;          % Number of sampling times
    inputs.exps.t_s{1}          = 0:5:duration;            % Times of samples

    inputs.exps.n_steps{1}      = sum(exps.t_con{1} < duration);
    inputs.exps.t_con{1}        = exps.t_con{1}(1:inputs.exps.n_steps{1}+1);
    inputs.exps.u{1}            = exps.u{1}(:,1:inputs.exps.n_steps{1});

    inputs.exps.exp_data{1}     = exps.exp_data{1}(1:inputs.exps.n_s{1},:);
    inputs.exps.error_data{1}   = exps.error_data{1}(1:inputs.exps.n_s{1},:);

    inputs.pathd.results_folder = results_folder;
    inputs.pathd.short_name     = short_name;
    inputs.pathd.runident       = strcat('pe-',int2str(i));

    % GLOBAL UNKNOWNS (SAME VALUE FOR ALL EXPERIMENTS)
    inputs.PEsol.id_global_theta=model.par_names;
    inputs.PEsol.global_theta_guess=transpose(global_theta_guess);
    inputs.PEsol.global_theta_max=global_theta_max;  % Maximum allowed values for the paramters
    inputs.PEsol.global_theta_min=global_theta_min;  % Minimum allowed values for the parameters

    % COST FUNCTION RELATED DATA
    inputs.PEsol.PEcost_type='llk';                       % 'lsq' (weighted least squares default) | 'llk' (log likelihood) | 'user_PEcost'
    inputs.PEsol.llk_type='homo_var';

    % SIMULATION
    inputs.ivpsol.ivpsolver='cvodes';
    inputs.ivpsol.senssolver='cvodes';
    inputs.ivpsol.rtol=1.0D-13;
    inputs.ivpsol.atol=1.0D-13;


    % OPTIMIZATION
    inputs.nlpsol.nlpsolver='eSS';
    inputs.nlpsol.eSS.maxeval = 2000;
    inputs.nlpsol.eSS.maxtime = 2000;
    inputs.nlpsol.eSS.local.solver = 'fmincon';  % nl2sol not yet installed on my mac
    inputs.nlpsol.eSS.local.finish = 'fmincon';  % nl2sol not yet installed on my mac
    inputs.rid.conf_ntrials=500;

    inputs.plotd.plotlevel='noplot';


    pe_start = now;
    results = AMIGO_PE(inputs);
    pe_inputs2{i} = inputs;
    pe_results2{i} = results;
    pe_end= now;

    % Write some results to the output file
    fid = fopen(resultFileName,'a');
    used_par_names = model.par_names;

    for j=1:size(used_par_names,1)
        fprintf(fid,'HOUR %d PARAM_FIT %s %f\n', i*5, used_par_names(j,:), results.fit.thetabest(j));
        if isfield(results.fit,'rel_conf_interval')
            fprintf(fid,'HOUR %d REL_CONF %s %f\n',  i*5, used_par_names(j,:), results.fit.rel_conf_interval(j));
        end
        if isfield(results.fit,'residuals')
           fprintf(fid,'HOUR %d RESIDUAL %s %f\n', i*5, used_par_names(j,:), results.fit.residuals{1}(j));
        end
        if isfield(results.fit,'rel_residuals')
            fprintf(fid,'HOUR %d REL_RESIDUAL %s %f\n', i*5, used_par_names(j,:), results.fit.rel_residuals{1}(j));
        end
    end
    % Time in seconds
    fprintf(fid,'HOUR %d PE_TIME %.1f\n',  i*5, (pe_end-pe_start)*24*60*60);
    fclose(fid);

    best_global_theta_log{i}=results.fit.thetabest;

end

true_param_values = model.par;

save([strcat(epccOutputResultFileNameBase),'.mat'], 'pe_inputs','pe_results','pe_results2','oed_results','exps','inputs','true_param_values','best_global_theta','best_global_theta_log');

out= 1;  
    

end