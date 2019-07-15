function [model] = ToggleSwitch_load_model_M1vs2_DVID()

% This script contains one of the model structures that we propose for the Toggle Switch
% Amongst the models being considered for model selection, this is Model
% M1.

% model.AMIGOjac = 0;                                                         % Compute Jacobian 0 = No, 1 = yes
model.input_model_type='charmodelC';                                        % Model introduction: 'charmodelC'|'c_model'|'charmodelM'|'matlabmodel'|'sbmlmodel'|'blackboxmodel'|'blackboxcost                             
model.n_st=8;                                                               % Number of states      
model.n_par=30;                                                             % Number of model parameters 
model.n_stimulus=2;                                                         % Number of inputs, stimuli or control variables   
model.stimulus_names=char('u_IPTG','u_aTc');                                % Name of stimuli or control variables
model.st_names=char('L_RFP','T_GFP','IPTGi','aTci',...
                    'L_RFP2','T_GFP2','IPTGi2','aTci2');      % Names of the states    
                
model.par_names=char('kL_p_m0','kL_p_m','theta_T','theta_aTc','n_aTc','n_T',...
                     'kT_p_m0','kT_p_m','theta_L','theta_IPTG','n_IPTG','n_L',...
                     'k_iptg',...
                     'k_aTc',...
                     'kL_p_m02','kL_p_m2','theta_T2','theta_aTc2','n_aTc2','n_T2',...
                     'kT_p_m02','kT_p_m2','theta_L2','theta_IPTG2','n_IPTG2','n_L2',...
                     'k_in_iptg','k_out_iptg',...
                     'k_in_aTc','k_out_aTc');                                  % Names of the parameters    
                 
model.eqns=...                                                              % Equations describing system dynamics. Time derivatives are regarded 'd'st_name''
               char('dL_RFP = 1/0.1386*(kL_p_m0 + (kL_p_m/(1+(T_GFP/theta_T*(1/(1+(aTci/theta_aTc)^n_aTc)))^n_T)))-0.0165*L_RFP',...
                    'dT_GFP = 1/0.1386*(kT_p_m0 + (kT_p_m/(1+(L_RFP/theta_L*(1/(1+(IPTGi/theta_IPTG)^n_IPTG)))^n_L)))-0.0165*T_GFP',...
                    'dIPTGi = k_iptg*(u_IPTG-IPTGi)-0.0165*IPTGi',...
                    'daTci = k_aTc*(u_aTc-aTci)-0.0165*aTci',...
                    'dL_RFP2 = 1/0.1386*(kL_p_m02 + (kL_p_m2/(1+(T_GFP2/theta_T2*(1/(1+(aTci2/theta_aTc2)^n_aTc2)))^n_T2)))-0.0165*L_RFP2',...
                    'dT_GFP2 = 1/0.1386*(kT_p_m02 + (kT_p_m2/(1+(L_RFP2/theta_L2*(1/(1+(IPTGi2/theta_IPTG2)^n_IPTG2)))^n_L2)))-0.0165*T_GFP2',...
                    'dIPTGi2 = k_iptg*(u_IPTG-IPTGi)-0.0165*IPTGi',...'%(k_in_IPTG*(u_IPTG-IPTGi2)+((k_in_IPTG*(u_IPTG-IPTGi2))^2)^0.5)/2-(k_out_IPTG*(IPTGi2-u_IPTG)+((k_out_IPTG*(IPTGi2-u_IPTG))^2)^0.5)/2',... 
                    'daTci2 = k_aTc*(u_aTc-aTci)-0.0165*aTci');%(k_in_aTc*(u_aTc-aTci2)+((k_in_aTc*(u_aTc-aTci2))^2)^0.5)/2-(k_out_aTc*(aTci2-u_aTc)+((k_out_aTc*(aTci2-u_aTc))^2)^0.5)/2');
                    

%==================
% PARAMETER VALUES
% =================

model.par=[0.05938212,6.27986,92.18167,0.9259555,0.4403347,3.295859,...
    0.2839974,4.49357,2.679409,0.01641903,0.7912365,0.7186982,...
    0.4464296,...
    0.01388927,...
    6.01684739e-02, 6.77117550e+00, 1.86956182e+01, 9.64252299e-01,7.72604972e-01, 2.66915511e+00,...
    1.23423001e-02, 5.13166491e+00, 4.69759475e+00, 6.05248321e-02, 1.02453826e+00, 5.21761158e-01,...
    1.70391261e-02, 1.83532118e-01, ...
    4.79131949e-02, 1.24263999e-02];

end               

        