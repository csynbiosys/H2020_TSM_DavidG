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
syms L0 L1 L2 Lac12 Lac12m G20 G21 G22 IPTGi Cit_molec

% unknown parameters
syms kLacI k2 kd km2 k1 km1 kLac12 kTP1 kcat Km kout kC lk        
 
% two outputs:
h    = [Cit_molec];   

% two inputs:
syms  IPTG
u    = [IPTG];

% states variables:
x    = [L0;L1;L2;Lac12;Lac12m;G20;G21;G22;IPTGi;Cit_molec];   

% parameters:
p    = [kLacI;k2;kd;km2;k1;km1;kLac12;kTP1;kcat;Km;kout;kC;lk];

% dynamic equations:
f    = [kLacI-(2*k2*IPTGi+kd)*L0+km2*L1-k1*(2*G20+G21)*L0+km1*(G21+2*G22);
        k2*(2*L0-L1)*IPTGi-(km2+kd)*L1+2*km2*L2;
        k2*L1*IPTGi-(2*km2+kd)*L2;
        kLac12-(kTP1+kd)*Lac12;
        kTP1*Lac12-kd*Lac12m;
        -2*k1*L0*G20+(km1+kd)*G21;
        2*k1*L0*G20+2*(km1+kd)*G22-(km1+k1*L0+kd)*G21;
        -2*(km1+kd)*G22+k1*L0*G21;
        (kcat*Lac12m*IPTG)/(Km+IPTG)-(kout*kd+2*k2*L0+k2*L1)*IPTGi+(kd+km2)*L1+2*(kd+km2)*L2;
        kC*G20+lk*kC*(G21+G22)-kd*Cit_molec];

% initial conditions: 
syms L00 L10 L20 Lac120 Lac12m0 G200 G210 G220 IPTGi0 Cit_molec0
ics = [];
% ics = [];

% which initial conditions are known:
% known_ics = [1 1 1 1 1 1]; 
known_ics = [0,0,0,0,0,0,0,0,0,1];

save('MIPr_model_DVID3','x','h','u','p','f','ics','known_ics');

