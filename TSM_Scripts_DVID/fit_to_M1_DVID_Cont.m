function [out] = fit_to_M1_DVID_Cont(epccOutputResultFileNameBase,epcc_exps,global_theta_guess)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit_to_ToggleSwitch - runs PE of the model structure we propose for the Toggle Switch to the calibration data published in Lugagne et al., 2017.
% PE is run, starting from 100 initial guesses for the parameters, on the whole dataset.
% An analytical solution for the steady state is used to define the initial
% conditions in a simulation of the system response to sustained inputs (48hrs)
% equal to the ON conditions. 
% Each vector of estimates is used to compute the SSE over the whole set. 
% Hence the estimate yielding the minimum SSE is selected. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Specify folder name and short_name
        results_folder = strcat('M1Fit',datestr(now,'yyyy-mm-dd-HHMMSS'));
        short_name     = strcat('M1F',int2str(epcc_exps));

        % Load experimental data. This data are derived, using the script
        % ExtractionAllExperimentalData/ExtractionStructure_AllData_final.m, starting from data used in the paper
        % Lugagne et al.
        % The data derives from 26 experiments, 6 calibration and 20 control
        % experiments.
        LugDat = load('D:\PhD\GitHub\H2020_Varun\AMIGO2_R2018b\Examples\TSM\model_selection\M2\AllDataLugagne_Final.mat');

        % Read the model into the model variable. Note that this model encodes
        % a path constraint to prevent negative solutions of ODEs. 
        % model=ToggleSwitch_load_model_M1(); 
        model=ToggleSwitch_load_model_M1_DVID_v2();

        % Initial guesses for theta 
        global_theta_min = [0.03*0.1, 10*0.1, 30*0.1, 0, 0, 0,...
            0.1*0.1, 1*0.1, 10*0.1, 0.01, 0, 0,...
            (5e-2)*0.1, (1e-1)*0.1]; 
        global_theta_max = [0.03*10, 10*10, 10*10, 100, 5, 5,...
            0.1*10, 1*10, 30*10, 1, 5, 5,...
            (5e-2)*10, (1e-1)*10];
        global_theta_guess = global_theta_guess;

        % Randomized vector of experiments
%         exps_indexall = 1:1:length(Data.expName);
        exps_indexall = [4,5,6,12,13,14,19,20,22,25];

        % Definition of the Training set. This has been obtained as
        % randperm(length(Data.expName)) and extracting the first 17 elements.
        % Note that here it is fixed to ensure that all global_theta_guess will
        % be trained on the same experiments.
%         exps_indexTraining = [22,6,3,16,11,7,17,14,8,5,21,25,2,19,15,1,23];
        exps_indexTraining = [4,6,13,19,22];

        % Definition of test set (remaining 9 experiments) 
%         exps_indexTest =  [2,4,18,24,13,9,20,10,12];
        exps_indexTest =  [5,12,14,20,25];

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Run  PE on the training set
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        best_global_theta = global_theta_guess; 

        % g_m, g_p, theta_T and theta_L excluded from identification
%         param_including_vector = [false,true,true,false,true,true,true,false,true,true,false,true,true,true,true,true,true,true];

        load(strcat(epccOutputResultFileNameBase,'PEMultiTest.mat'));

        %%    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Run simulation on the test set
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear inputs;
        clear exps;
        inputs.model = model;
        inputs.model.par = best_global_theta;

        inputs.pathd.results_folder = results_folder;                        
        inputs.pathd.short_name     = short_name;
        inputs.pathd.runident       = 'initial_setup2';

        AMIGO_Prep(inputs);

        inputs.pathd.results_folder = results_folder;                        
        inputs.pathd.short_name     = short_name;
        inputs.pathd.runident       = strcat('simTest-',int2str(epcc_exps));


        % Compute the steady state considering the initial theta guess, u_IPTG and
        % u_aTc
        y0_test = zeros(length(exps_indexTraining),model.n_st);

        for iexp=1:length(exps_indexTraining)
            inexp = exps_indexTraining(iexp);
            y0_test(iexp,:) = M1_Compute_SteadyState_OverNight_DVID_v2(epcc_exps,inputs,1,global_theta_guess,...
                [LugDat.Data.exp_data{inexp}(1,1), LugDat.Data.exp_data{inexp}(2,1)],...
                [LugDat.Data.Initial_IPTG{inexp} LugDat.Data.Initial_aTc{inexp}]+1e-7);
        end

        exps.n_exp = length(exps_indexTest);

        for iexp=1:length(exps_indexTest)

            inexp = exps_indexTest(iexp);
            duration = round(LugDat.Data.t_samples{inexp}(1,end));                                                
            inducer = LugDat.Data.input{inexp}+1e-7;
            
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
            exps.exp_y0{iexp} = y0_test(iexp,:);

        end


        inputs.exps  = exps;

        inputs.pathd.results_folder = results_folder;
        inputs.pathd.short_name     = short_name;
        inputs.pathd.runident       ='-sim';
        inputs.plotd.plotlevel='min';

        sim_results = AMIGO_SObs(inputs);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%Compute SSE on the test set
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SSE = zeros(length(exps_indexTest),2);
        for iexp=1:length(exps_indexTest)
            exp_indexData = exps_indexTest(iexp);
            SSE_RFP(iexp) = sum((LugDat.Data.exp_data{1,exp_indexData}(1,:)-sim_results.sim.sim_data{1,iexp}(:,1)').^2);
            SSE_GFP(iexp) = sum((LugDat.Data.exp_data{1,exp_indexData}(2,:)-sim_results.sim.sim_data{1,iexp}(:,2)').^2);
            SSE(iexp,:) = [SSE_RFP(iexp),SSE_GFP(iexp)];
        end


        sim_inputs = inputs;
        sim_exps = exps;
        save(strcat(epccOutputResultFileNameBase,'-sim','.mat'),'sim_results','sim_inputs','sim_exps','best_global_theta','SSE');

out = 1;
end