function [] = Run_fit_to_M1_DVID( resultBase, numExperiments )
% This function runs the script for the identification of the model structure
% we propose to describe the Toggle switch to the experimental data published in Lugagne et al., 2017.
% The dataset includes 26 experiments (6 calibration+20 control). Inconsistencies noted in the original repository of
% the paper have been addressed through private communication with the
% authors and corrected accordingly. Parameter estimation is performed
% using crossvalidation. To this end we split the data between training set
% (17 randomised experiments)and test (9 experiments) set. SSE is computed
% over the test set to select the best parameter estimates. 
% It takes two inputs:
% - a string with the identifier of the output files
% - the number of iterations of parameter estimation. 
% In the paper we used numExperiments = 100.

% cd ('../../../../../');
% AMIGO_Startup();
% cd('Examples\ToggleSwitch_FM_LB_VBK\model_calibration\M1\Parameter_estimation');

% Create a matrix of initial guesses for the parameters, having as many
% rows as the number of PE iterations (numExperiments)
% Each vector is passed as input to the computing function

% theta_min = [0.1386,5e-06,1e-3,40,1,2,2,0.0165,1e-05,0.001,0.0215,0.0001,2,2,10,10,0.001,0.001]; % verify Theta_T is correct 
% theta_max = [0.1386,10,1000,40,100,4,4,0.0165,30,1000,0.0215,0.01,4,4,1000,1000,0.1,0.1];
% 
% M_norm = lhsdesign(numExperiments,length(theta_min));
% M = zeros(size(M_norm));
% for c=1:size(M_norm,2)
%     for r=1:size(M_norm,1)
%         M(r,c) = 10^(M_norm(r,c)*(log10(theta_max(1,c))-log10(theta_min(1,c)))+log10(theta_min(1,c))); % log exploration
%     end
% end 
% 
% % %check the location of the parameters that are fixed
%  ParFull = M; %[0.1386*ones(size(M,1),1) M(:,1:6) 0.0165*ones(size(M,1),1) M(:,7:end)];
%  
% save('MatrixParameters.mat','ParFull');


% global_theta_guess = [0.05938212,6.27986,92.18167,0.9259555,0.4403347,3.295859,0.2839974,4.49357,2.679409,...
%     0.01641903,0.7912365,0.7186982,0.4464296,0.01388927];

% ParFull = load('draws_ALL_Long_Model3.stan.mat');
% c = ParFull.x(:,3:14);
% c = [c ParFull.x(:,1:2)];

ParFull = load('MatrixParameters.mat');

parfor epcc_exps=1:numExperiments
        try
            global_theta_guess = ParFull.ParFull(epcc_exps,:);
            epccOutputResultFileNameBase = [resultBase,'-',num2str(epcc_exps)];
            [out] = fit_to_M1_DVID(epccOutputResultFileNameBase,epcc_exps,global_theta_guess);
        catch err
            %open file
            errorFile = [resultBase,'-',num2str(epcc_exps),'.errorLog'];
            fid = fopen(errorFile,'a+');
            fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'));
            % close file
            fclose(fid);
        end
end

