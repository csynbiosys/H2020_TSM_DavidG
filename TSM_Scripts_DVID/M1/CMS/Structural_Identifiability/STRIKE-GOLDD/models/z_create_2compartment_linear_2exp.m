% Markov model, 2 states, linear chain, with u=constant and 2 experiments:
clear all;

% 2 -> 4 states
syms x1Exp1 x2Exp1 x1Exp2 x2Exp2
x = [x1Exp1;x2Exp1;x1Exp2;x2Exp2];

% 1 -> 2 outputs
h = [x1Exp1;x1Exp2];

% 1 -> 2 inputs
syms u1Exp1 u1Exp2;
u = [u1Exp1;u1Exp2];
% u = [];

% 4 unknown parameters 
syms x3 x4 x5 x6
p =[x3; x4; x5; x6];

% dynamic equations
f = [-(x3+x4)*x1Exp1+x5*x2Exp1+x6*u1Exp1;  % Exp1
     x4*x1Exp1-x5*x2Exp1;                  % Exp1
     -(x3+x4)*x1Exp2+x5*x2Exp2+x6*u1Exp2;  % Exp2
     x4*x1Exp2-x5*x2Exp2];                 % Exp2;

% initial conditions
ics  = [0.5,0.5,0.5,0.5]; 
known_ics = [1,0,1,0];

save('two_compartment_linear_2Exp','x','p','h','f','u','ics','known_ics');

