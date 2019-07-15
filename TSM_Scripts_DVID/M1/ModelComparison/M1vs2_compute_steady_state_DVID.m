function [ res ] = M1vs2_compute_steady_state_DVID(theta, InitialStates_AU, initial_u)
% ToggleSwitch_M1_Compute_SteadyState computes the steady state of the MToggleSwitch model for the given values of
% theta and the inputs u_IPTG and u_aTc.

% Model 1
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
k_iptg = theta(13);
k_aTc = theta(14);

% Model 2
kL_p_m02 = theta(15);
kL_p_m2 = theta(16);
theta_T2 = theta(17);
theta_aTc2 = theta(18);
n_aTc2 = theta(19);
n_T2 = theta(20);
kT_p_m02 = theta(21);
kT_p_m2 = theta(22);
theta_L2 = theta(23);
theta_IPTG2 = theta(24);
n_IPTG2 = theta(25);
n_L2 = theta(26);
k_in_iptg = theta(27);
k_out_iptg = theta(28);
k_in_aTc = theta(29);
k_out_aTc = theta(30);



L_RFP = InitialStates_AU(1);
T_GFP = InitialStates_AU(2);
L_RFP2 = InitialStates_AU(1);
T_GFP2 = InitialStates_AU(2);

u_IPTG = initial_u(1);
u_aTc = initial_u(2);

%% Steady state equation

aTci = k_aTc*u_aTc/(0.0165+k_aTc);
IPTGi = k_iptg*u_IPTG/(0.0165+k_iptg);

L_molec = (1/(0.0165*0.1386))*(kL_p_m0 + (kL_p_m/(1+(T_GFP/theta_T*(1/(1+(aTci/theta_aTc)^n_aTc)))^n_T)));

T_molec = (1/(0.0165*0.1386))*(kT_p_m0 + (kT_p_m/(1+(L_RFP/theta_L*(1/(1+(IPTGi/theta_IPTG)^n_IPTG)))^n_L)));

aTci2 = u_aTc;
IPTGi2 = u_IPTG;

L_molec2 = (1/(0.0165*0.1386))*(kL_p_m02 + (kL_p_m2/(1+(T_GFP2/theta_T2*(1/(1+(aTci2/theta_aTc2)^n_aTc2)))^n_T2)));

T_molec2 = (1/(0.0165*0.1386))*(kT_p_m02 + (kT_p_m2/(1+(L_RFP2/theta_L2*(1/(1+(IPTGi2/theta_IPTG2)^n_IPTG2)))^n_L2)));

res = [L_molec T_molec IPTGi aTci L_molec2 T_molec2 IPTGi2 aTci2];

end
