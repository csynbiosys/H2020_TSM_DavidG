%--------------------------------------------------------------------------
% File that creates the Toggle Switch model. 
% It stores it in a mat-file named TSM.mat.
% 
%--------------------------------------------------------------------------
% Re-formulation of the toggle switch model described in  
% J.-B. Lugagne, S. S. Carrillo, M. Kirch, A. Köhler, G. Batt, and P. Hersen,
% “Balancing a genetic toggle switch by real-time feedback control and periodic
% forcing,” Nature Communications, vol. 8, no. 1, p. 1671, Nov. 2017.
%--------------------------------------------------------------------------
clear all;

% states
syms L_RFP T_GFP IPTGi aTci

% unknown parameters
syms kL_p_m0 kL_p_m theta_T theta_aTc n_aTc n_T ...
         kT_p_m0 kT_p_m theta_L theta_IPTG n_IPTG n_L...
         k_iptg k_aTc         
 
% two outputs:
h    = [L_RFP;T_GFP];   

% two inputs:
syms  u_IPTG u_aTc
u    = [u_IPTG; u_aTc];

% states variables:
x    = [L_RFP;T_GFP;IPTGi;aTci];   

% parameters:
p    = [kL_p_m0;kL_p_m;theta_T;theta_aTc;n_aTc;n_T;...
         kT_p_m0;kT_p_m;theta_L;theta_IPTG;n_IPTG;n_L;...
         k_iptg;k_aTc];

% dynamic equations:
f    = [1/0.1386*(kL_p_m0 + (kL_p_m/(1+(T_GFP/theta_T*(1/(1+(aTci/theta_aTc)^n_aTc)))^n_T)))-0.0165*L_RFP;
        1/0.1386*(kT_p_m0 + (kT_p_m/(1+(L_RFP/theta_L*(1/(1+(IPTGi/theta_IPTG)^n_IPTG)))^n_L)))-0.0165*T_GFP;
        k_iptg*(u_IPTG-IPTGi)-0.0165*IPTGi;
        k_aTc*(u_aTc-aTci)-0.0165*aTci];
%         k_aTc*(u_aTc-aTci)-0.0165*aTci];  

% initial conditions: 
syms L_RFP0 T_GFP0 IPTGi0 aTci0
% ics = [L_RFP0,T_GFP0,IPTGi0,aTci0];
ics = [];

% which initial conditions are known:
% known_ics = [1 1 1 1 1 1]; 
known_ics = [0,0,0,0];

save('TSM_model1_DVID','x','h','u','p','f','ics','known_ics');

