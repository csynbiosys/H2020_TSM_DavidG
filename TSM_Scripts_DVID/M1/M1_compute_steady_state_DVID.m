function [ res ] = M1_compute_steady_state_DVID(theta, InitialStates_AU, initial_u)
% ToggleSwitch_M1_Compute_SteadyState computes the steady state of the MToggleSwitch model for the given values of
% theta and the inputs u_IPTG and u_aTc.

kL_p_m0 = theta(1);
kL_p_m = theta(2);
theta_T = theta(3);
theta_aTc = theta(4);
n_aTc = theta(5);
n_T = theta(6);
kT_p_m0 = theta(7);
kT_p_m = theta(8);
theta_L = theta(9);
theta_IPTG = theta(10);
n_IPTG = theta(11);
n_L = theta(12);
sc_T_molec = theta(13);
sc_L_molec = theta(14);
k_iptg = theta(15);
k_aTc = theta(16);

L_AU = InitialStates_AU(1);
T_AU = InitialStates_AU(2);

u_IPTG = initial_u(1);
u_aTc = initial_u(2);

%% Steady state equation

L_molec = L_AU/sc_L_molec;
T_molec = T_AU/sc_T_molec;

aTci = k_aTc*u_aTc/(0.0165+k_aTc);
IPTGi = k_iptg*u_IPTG/(0.0165+k_iptg);

L_molec = (1/(0.0165*0.1386))*(kL_p_m0 + (kL_p_m/(1+(T_molec/theta_T*(1/(1+(aTci/theta_aTc)^n_aTc)))^n_T)));

T_molec = (1/(0.0165*0.1386))*(kT_p_m0 + (kT_p_m/(1+(L_molec/theta_L*(1/(1+(IPTGi/theta_IPTG)^n_IPTG)))^n_L)));

L_AU = sc_L_molec*L_molec;
T_AU = sc_T_molec*T_molec;


res = [L_molec T_molec T_AU L_AU IPTGi aTci];

end

