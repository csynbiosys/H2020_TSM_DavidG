function [ initial_values ] = M1_compute_steady_state(par, initial_y, initial_u)
%ToggleSwitch_Compute_SteadyState Calculates the steady state for par,
%u_IPTG and u_aTc.
%   Computes the steady state of the MToggleSwitch model for the given values of
%   par and the inputs u_IPTG and u_aTc.

g_m = par(1);
kL_p_m0 = par(2);
kL_p_m = par(3);
par_T = par(4);
par_aTc = par(5);
n_aTc = par(6);
n_T = par(7);
g_p = par(8);
kT_p_m0 = par(9);
kT_p_m = par(10);
par_L = par(11);
par_IPTG = par(12);
n_IPTG = par(13);
n_L = par(14);
sc_T_molec = par(15);
sc_L_molec = par(16);
k_iptg = par(17);
k_aTc = par(18);

T_AU = initial_y(3); % Get from initial conditions
L_AU = initial_y(4); % EBC changed this- prvious version TAU y(4) LAU y(3) --> confirm with Varun & Lucia



% Inputs
u_IPTG=initial_u(1); 
u_aTc=initial_u(2);
 
%% Steady state equation

L_molec1 = L_AU/sc_L_molec;
T_molec1 = T_AU/sc_T_molec;

aTci = k_aTc*u_aTc/(g_p+k_aTc);
IPTGi = k_iptg*u_IPTG/(g_p+k_iptg);

L_molec = (1/(g_p*g_m))*(kL_p_m0 + (kL_p_m/(1+(T_molec1/par_T*(1/(1+(aTci/par_aTc)^n_aTc)))^n_T)));

T_molec = (1/(g_p*g_m))*(kT_p_m0 + (kT_p_m/(1+(L_molec1/par_L*(1/(1+(IPTGi/par_IPTG)^n_IPTG)))^n_L)));

L_AU = sc_L_molec*L_molec; 
T_AU = sc_T_molec*T_molec; 


initial_values = [L_molec T_molec T_AU L_AU IPTGi aTci 0];  


end

