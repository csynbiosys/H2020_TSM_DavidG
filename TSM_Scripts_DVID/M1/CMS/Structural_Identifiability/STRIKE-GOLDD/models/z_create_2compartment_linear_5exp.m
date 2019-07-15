% Markov model, 2 states, linear chain, with u=constant and 2 experiments:
clear all;

% 2 -> 10 states
syms x1Exp1 x2Exp1 x1Exp2 x2Exp2 x1Exp3 x2Exp3 x1Exp4 x2Exp4 x1Exp5 x2Exp5
x = [x1Exp1 x2Exp1 x1Exp2 x2Exp2 x1Exp3 x2Exp3 x1Exp4 x2Exp4 x1Exp5 x2Exp5].';

% 1 -> 4 outputs
h = [x1Exp1 x1Exp2 x1Exp3 x1Exp4 x1Exp5].';

% 1 -> 4 inputs
syms u1Exp1 u1Exp2 u1Exp3 u1Exp4 u1Exp5
u = [u1Exp1 u1Exp2 u1Exp3 u1Exp4 u1Exp5].';

% 4 unknown parameters 
syms x3 x4 x5 x6
p =[x3; x4; x5; x6];

% dynamic equations
f = [-(x3+x4)*x1Exp1+x5*x2Exp1+x6*u1Exp1;  % Exp1
     x4*x1Exp1-x5*x2Exp1;                  % Exp1
     -(x3+x4)*x1Exp2+x5*x2Exp2+x6*u1Exp2;  % Exp2
     x4*x1Exp2-x5*x2Exp2;                  % Exp2
     -(x3+x4)*x1Exp3+x5*x2Exp3+x6*u1Exp3;  % Exp3
     x4*x1Exp3-x5*x2Exp3;                  % Exp3
     -(x3+x4)*x1Exp4+x5*x2Exp4+x6*u1Exp4;  % Exp4
     x4*x1Exp4-x5*x2Exp4;                  % Exp4
     -(x3+x4)*x1Exp5+x5*x2Exp5+x6*u1Exp5;  % Exp5
     x4*x1Exp5-x5*x2Exp5];                 % Exp5 
% initial conditions
ics  = []; 
known_ics = [1,0,1,0,1,0,1,0,1,0];

save('two_compartment_linear_5Exp','x','p','h','f','u','ics','known_ics');

