function [] = run_in_silico_experiment_parfor_IntuitionDriven_M1_DVID(resultBase, numExperiments )

% cd ('../../../../../');
% AMIGO_Startup();
% 
% cd ('/Examples/TSM/model_calibration/M1/CIC');

% Upper and lower boundaries selected for the parameters (theta refers to the parameter vector)
% theta_min=[0.1386,5e-06,1e-3,7000,1,2,2,0.0165,1e-05,0.001,700,3.07e-4,2,2,1.32,1,0.001,0.001]; % verify Theta_T is correct 
% theta_max=[0.1386,10,1000,7000,100,4,4,0.0165,30,1000,700,3.07e-2,4,4,1320,1000,0.1,0.1]; % Maximum allowed values for the parameters
          
% Create a matrix of initial guesses for the vector of parameter estimates.
% The matrix has as many rows as the number of PE iterations
% (numExperiments).
% The vectors are obtained as latin hypercube samples within the boundaries selected above.
% Each vector is passed as input to the computing function

% M_norm = lhsdesign(numExperiments,length(theta_min));
% M = zeros(size(M_norm));
% for c=1:size(M_norm,2)
%     for r=1:size(M_norm,1)
%         M(r,c) = 10^(M_norm(r,c)*(log10(theta_max(1,c))-log10(theta_min(1,c)))+log10(theta_min(1,c))); % log exploration
%     end
% end 
% 
% ParFull = M; % in this case I am fitting all the values
% save('MatrixParameters_InputComparison.mat','ParFull');

% For the sake of a fair comparison among input classes, all simulations
% used the same matrix of initial theta guesses
% load('MatrixParameters_InputComparison.mat');       
% load('AllDataLugagne_Final_21exp.mat');

parfor epcc_exps=1:numExperiments
        try
            global_theta_guess = [0.05938212,6.27986,92.18167,0.9259555,0.4403347,3.295859,0.2839974,4.49357,2.679409,0.01641903,0.7912365,0.7186982,505,505,0.4464296,0.01388927];%ParFull(epcc_exps,:);
            epccOutputResultFileNameBase = [resultBase,'-',num2str(epcc_exps)];
            
            % Uncomment the desired function for parameter estimation
            %[out] = fit_to_Toggle_switch_Step(epccOutputResultFileNameBase,epcc_exps,global_theta_guess); % How uses the homovar data (5%)
            %[out] = fit_to_Toggle_switch_Pulse(epccOutputResultFileNameBase,epcc_exps,global_theta_guess); % How uses the homovar data (5%)
            [out] = fit_to_M1_Random_DVID(epccOutputResultFileNameBase,epcc_exps,global_theta_guess); % How uses the homovar data (5%)

        catch err
            %open file
            errorFile = [resultBase,'-',num2str(epcc_exps),'.errorLog'];
            fid = fopen(errorFile,'a+');
            fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'));
            % close file
            fclose(fid);
        end
end
