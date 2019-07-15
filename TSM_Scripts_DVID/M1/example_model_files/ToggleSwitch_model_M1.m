function [model] = ToggleSwitch_model_M1()

% This script contains one of the model structures that we propose for the Toggle Switch
% Amongst the models being considered for model selection, this is Model
% M1.

model.AMIGOjac = 0;                                                         % Compute Jacobian 0 = No, 1 = yes
model.input_model_type='charmodelC';                                        % Model introduction: 'charmodelC'|'c_model'|'charmodelM'|'matlabmodel'|'sbmlmodel'|'blackboxmodel'|'blackboxcost                             
model.n_st=7;                                                               % Number of states      
model.n_par=18;                                                             % Number of model parameters 
model.n_stimulus=2;                                                         % Number of inputs, stimuli or control variables   
model.stimulus_names=char('u_IPTG','u_aTc');                                % Name of stimuli or control variables
model.st_names=char('L_molec','T_molec','T_AU','L_AU','IPTGi','aTci','sviol');     % Names of the states                                              
model.par_names=char('g_m','kL_p_m0','kL_p_m','theta_T','theta_aTc','n_aTc','n_T','g_p',...
                     'kT_p_m0','kT_p_m','theta_L','theta_IPTG','n_IPTG','n_L',...
                     'sc_T_molec',...
                     'sc_L_molec',...
                     'k_iptg',...
                     'k_aTc');                                 % Names of the parameters    
                 
model.eqns=...                                                              % Equations describing system dynamics. Time derivatives are regarded 'd'st_name''
               char('dL_molec = 1/g_m*(kL_p_m0 + (kL_p_m/(1+(T_molec/theta_T*(1/(1+(aTci/theta_aTc)^n_aTc)))^n_T)))-g_p*L_molec',...
                    'dT_molec = 1/g_m*(kT_p_m0 + (kT_p_m/(1+(L_molec/theta_L*(1/(1+(IPTGi/theta_IPTG)^n_IPTG)))^n_L)))-g_p*T_molec',...
                    'dT_AU = sc_T_molec*dT_molec',...
                    'dL_AU = sc_L_molec*dL_molec',...
                    'dIPTGi = k_iptg*(u_IPTG-IPTGi)-g_p*IPTGi',...
                    'daTci = k_aTc*(u_aTc-aTci)-g_p*aTci',...
                    'dsviol = (fmax(-L_molec,0))^2+ (fmax(-T_molec,0))^2+(fmax(-IPTGi,0))^2+ (fmax(-aTci,0))^2'); % single constraint to be handled by the optimiser. If violated (e.g. one of the states becomes negative) the constraint becomes greater than 0

%==================
% PARAMETER VALUES
% =================

        
g_m=0.1386;       kL_p_m0=4.2083e-1;
kL_p_m=2.5235e+1;    theta_T=40;
theta_aTc=3.7453e-01; n_aTc=1.0415e+00;
n_T=1.6342e+00;   g_p=0.01650;

kT_p_m0=1.1161e+02;   kT_p_m=1.0645e+03;
theta_L=0.0222;   theta_IPTG=6.4438e-01;
n_IPTG=3.2915e+00;    n_L=1.1361e+00;
sc_T_molec=2.8962e-03;    sc_L_molec=3.5841e-01;
k_iptg=7.1516e-02;    k_aTc= 7.8176e-03;



model.par=[g_m,kL_p_m0,kL_p_m,theta_T,theta_aTc,n_aTc,n_T,g_p,...
                     kT_p_m0,kT_p_m,theta_L,theta_IPTG,n_IPTG,n_L,...
                     sc_T_molec,...
                     sc_L_molec,...
                     k_iptg,...
                     k_aTc];

end                                 