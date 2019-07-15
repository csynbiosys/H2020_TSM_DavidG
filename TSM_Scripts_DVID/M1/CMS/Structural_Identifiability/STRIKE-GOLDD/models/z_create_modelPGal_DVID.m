%--------------------------------------------------------------------------
% File that creates the MIPr model. 
% It stores it in a mat-file named MIPr.mat.
% 
%--------------------------------------------------------------------------
% Re-formulation of the toggle switch model described in  
% J.-B. Lugagne, S. S. Carrillo, M. Kirch, A. Köhler, G. Batt, and P. Hersen,
% “Balancing a genetic toggle switch by real-time feedback control and periodic
% forcing,” Nature Communications, vol. 8, no. 1, p. 1671, Nov. 2017.
%--------------------------------------------------------------------------
clear all;

% states
syms Cit_mrna Cit_foldedP Cit_fluo

% unknown parameters
syms alpha1 Vm1 h1 Km1 d1 alpha2 d2 Kf         
 
% two outputs:
h    = [Cit_fluo];   

% two inputs:
syms  IPTG
u    = [IPTG];

% states variables:
x    = [Cit_mrna;Cit_foldedP;Cit_fluo];   

% parameters:
p    = [alpha1;Vm1;h1;Km1;d1;alpha2;d2;Kf];

% dynamic equations:
f    = [alpha1+(Vm1*(IPTG^h1/(Km1^h1+IPTG^h1)))-d1*Cit_mrna;
        alpha2*Cit_mrna-(d2+Kf)*Cit_foldedP;
        Kf*Cit_foldedP-d2*Cit_fluo];

% initial conditions: 
syms Cit_mrna0 Cit_foldedP0 Cit_fluo0
% ics = [0 0 0];
ics = [];

% which initial conditions are known:
% known_ics = [1 1 1 1 1 1]; 
known_ics = [0 0 1];

save('MIPr_model_DVID','x','h','u','p','f','ics','known_ics');

