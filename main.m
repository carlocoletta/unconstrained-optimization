clc
clear all

%The main is structured in different sections to be run modifing the
%arguments to adapt at every problem the requests

%Initial settings
random_seed = min(294507, 344810);
rng(random_seed);
dimensions_newton = [1e3, 1e4, 1e5];
dimensions_nelder = [10, 25, 50];

%Nelder parameters: we pass them as a vector
rho = 1.2;
chi = 2.2;
gamma = 0.5;
sigma = 0.5;
kmax = 5000;
tol = 1e-4;

parameters_nelder = [rho, chi, gamma, sigma, kmax, tol];


%Newton parameters: we pass them as a vector
rho = 0.5;
c1 = 1e-4;
btmax = 50;
tolgrad = 1e-4; 
tau_kmax = 100;
kmax = 500;

parameters_newton = [rho, c1, btmax, tolgrad, tau_kmax, kmax];


%% 
%FUNCTION PLOTTING
f = problema1();
generalFunctionPlot(f, [0,2], [0,2], 100, "Extended Rosenbrock");

%% 
%ROSENBROCK FUNCTION
rosenbrockInitialCheck();

%%
%NELDER MEAD
wrapperNelder(dimensions_nelder, 1, parameters_nelder);

%%
%NEWTON WITH EXACT DERIVATIVES
wrapperNewtonExact(dimensions_newton, 1, parameters_newton);

%% 
%NEWTON WITH FINITE DIFFERENCES AND FIXED STEP (FLAG 0)
for i=1:3
    wrapperNewtonFinDiff(dimensions_newton(i), 1, parameters_newton, 0);
end


%%
%NEWTON WITH FINITE DIFFERENCES AND VARIABLE STEP (FLAG 1)
for i=1:3
    wrapperNewtonFinDiff(dimensions_newton(i), 1, parameters_newton, 1);
end





