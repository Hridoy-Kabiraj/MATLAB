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
N0 = [1/lambda_Xe138; 0];  % Use column vector

%% Time span (minutes)
tSpan = [0 400];

%% Solve with ode45
[t_ode45, N_ode45] = ode45(decayChain, tSpan, N0);
% Compute activities
A_Xe_ode45 = lambda_Xe138 * N_ode45(:,1);
A_Cs_ode45 = lambda_Cs138 * N_ode45(:,2);

%% Solve with rk4_solver
h = 0.5; % Define step size for RK4
[t_rk4, N_rk4] = rk4_solver(decayChain, tSpan, N0, h);
% Compute activities
A_Xe_rk4 = lambda_Xe138 * N_rk4(:,1);
A_Cs_rk4 = lambda_Cs138 * N_rk4(:,2);

figure;
% ode45 results 
plot(t_ode45, A_Xe_ode45, 'm-', 'LineWidth', 2, 'DisplayName', '^{138}Xe (ode45)'); hold on;
plot(t_ode45, A_Cs_ode45, 'b-', 'LineWidth', 2, 'DisplayName', '^{138}Cs (ode45)');

% rk4 results 
plot(t_rk4, A_Xe_rk4, 'm--', 'LineWidth', 1.5, 'DisplayName', '^{138}Xe (RK4)');
plot(t_rk4, A_Cs_rk4, 'b--', 'LineWidth', 1.5, 'DisplayName', '^{138}Cs (RK4)');
hold off;

xlabel('Time (minutes)');
ylabel('Activity (decays/min)');
legend('Location', 'northeast');
title('RK4 vs ode45 Comparison: ^{138}Xe → ^{138}Cs');
subtitle('No Equilibrium');
grid on;