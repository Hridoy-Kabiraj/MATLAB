clc; clearvars; close all;

%% Decay constants (half-lives in minutes)
lambda_Xe138 = log(2) / 14.08;   % Xe-138
lambda_Cs138 = log(2) / 33.41;   % Cs-138

%% ODE system for the decay chain
decayChain = @(t, N) [
    -lambda_Xe138 * N(1);                % Xe-138
     lambda_Xe138 * N(1) - lambda_Cs138 * N(2)  % Cs-138
];

%% Initial conditions
N0 = [1/lambda_Xe138, 0];  % Initial number of nuclei

%% Time span (minutes)
tSpan = [0 400];

%% Solve the ODE system
[t, N] = ode45(decayChain, tSpan, N0);

%% Plot Activity vs Time
figure;
plot(t, lambda_Xe138*N(:,1), 'm-', ...
     t, lambda_Cs138*N(:,2), 'b-', 'LineWidth', 1.5);

xlabel('Time (minutes)');
ylabel('Activity (decays/min)');
legend('^{138}Xe (Parent)', '^{138}Cs (Daughter)', 'Location', 'northeast');
title('Radioactive Decay of ^{138}Xe → ^{138}Cs');
subtitle('No Equilibrium');
grid on;
