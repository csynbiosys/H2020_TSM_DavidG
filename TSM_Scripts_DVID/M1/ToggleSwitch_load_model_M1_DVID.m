function [model] = ToggleSwitch_load_model_M1_DVID()

% This script contains one of the model structures that we propose for the Toggle Switch
% Amongst the models being considered for model selection, this is Model
% M1.

model.AMIGOjac = 0;                                                         % Compute Jacobian 0 = No, 1 = yes
model.input_model_type='charmodelC';                                        % Model introduction: 'charmodelC'|'c_model'|'charmodelM'|'matlabmodel'|'sbmlmodel'|'blackboxmodel'|'blackboxcost                             
model.n_st=6;                                                               % Number of states      
model.n_par=16;                                                             % Number of model parameters 
model.n_stimulus=2;                                                         % Number of inputs, stimuli or control variables   
model.stimulus_names=char('u_IPTG','u_aTc');                                % Name of stimuli or control variables
model.st_names=char('L_molec','T_molec','T_AU','L_AU','IPTGi','aTci');      % Names of the states                                              
model.par_names=char('kL_p_m0','kL_p_m','theta_T','theta_aTc','n_aTc','n_T',...
                     'kT_p_m0','kT_p_m','theta_L','theta_IPTG','n_IPTG','n_L',...
                     'sc_T_molec',...
                     'sc_L_molec',...
                     'k_iptg',...
                     'k_aTc');                                  % Names of the parameters    
                 
model.eqns=...                                                              % Equations describing system dynamics. Time derivatives are regarded 'd'st_name''
               char('dL_molec = 1/0.1386*(kL_p_m0 + (kL_p_m/(1+(T_molec/theta_T*(1/(1+(aTci/theta_aTc)^n_aTc)))^n_T)))-0.0165*L_molec',...
                    'dT_molec = 1/0.1386*(kT_p_m0 + (kT_p_m/(1+(L_molec/theta_L*(1/(1+(IPTGi/theta_IPTG)^n_IPTG)))^n_L)))-0.0165*T_molec',...
                    'dT_AU = sc_T_molec*dT_molec',...
                    'dL_AU = sc_L_molec*dL_molec',...
                    'dIPTGi = k_iptg*(u_IPTG-IPTGi)-0.0165*IPTGi',...
                    'daTci = k_aTc*(u_aTc-aTci)-0.0165*aTci');
                    

%==================
% PARAMETER VALUES
% =================

% dummy parameter values set to the mid point of the range of each individual parameter value
% model.par=[0.1386,5.0000025,500.0005,40,0.213428215,3,3,0.01650,15.000005,500.0005,0.0222,0.00505,3,3,505,505,0.0505,0.0505]; % Parameter values set to the average of global_theta_guess min, max
                 
% Parameters from Eva's report (9/11/2018 16:13) on PE carried out using experiments 2, 19
% and 15. 

model.par=[0.05938212,6.27986,92.18167,0.9259555,0.4403347,3.295859,0.2839974,4.49357,2.679409,...
    0.01641903,0.7912365,0.7186982,1,1,0.4464296,0.01388927];

end                                 