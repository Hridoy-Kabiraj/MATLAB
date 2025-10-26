clc; clearvars; close all;

%% Decay constants (half-lives in hours)
lambda_Mo99  = log(2)/66;        % 99Mo
lambda_Tc99m = log(2)/6;         % 99mTc
lambda_Tc99  = log(2)/210000;    % 99Tc

%% ODE system for the decay chain
decayChain = @(t, N) [
    -lambda_Mo99 * N(1);                           % Mo-99
    0.88 * lambda_Mo99 * N(1) - lambda_Tc99m * N(2); % Tc-99m
    0.12 * lambda_Mo99 * N(1) - lambda_Tc99 * N(3)   % Tc-99
];

%% Initial conditions
N0 = [1/lambda_Mo99; 0; 0];   % Use column vector

%% Time span
tSpan = [0, 100];  % hours

%% Solve with ode45
% Use a finer tSpan for ode45 output if needed, but [t0, tf] is fine for solving
[t_ode45, N_ode45] = ode45(decayChain, tSpan, N0);

%% Solve with rk4_solver
h = 0.1; % Define step size for RK4
[t_rk4, N_rk4] = rk4_solver(decayChain, tSpan, N0, h);

figure;
% ode45 results
plot(t_ode45, lambda_Mo99*N_ode45(:,1), 'm-', 'LineWidth', 2, 'DisplayName', '^{99}Mo (ode45)'); hold on;
plot(t_ode45, lambda_Tc99m*N_ode45(:,2), 'c-', 'LineWidth', 2, 'DisplayName', '^{99m}Tc (ode45)');
plot(t_ode45, lambda_Tc99*N_ode45(:,3), 'b-', 'LineWidth', 2, 'DisplayName', '^{99}Tc (ode45)');

% rk4 results
plot(t_rk4, lambda_Mo99*N_rk4(:,1), 'm--', 'LineWidth', 1.5, 'DisplayName', '^{99}Mo (RK4)');
plot(t_rk4, lambda_Tc99m*N_rk4(:,2), 'c--', 'LineWidth', 1.5, 'DisplayName', '^{99m}Tc (RK4)');
plot(t_rk4, lambda_Tc99*N_rk4(:,3), 'b--', 'LineWidth', 1.5, 'DisplayName', '^{99}Tc (RK4)');
hold off;

xlabel('Time (hours)');
ylabel('Activity (decays/hour)');
legend('Location', 'best');
title('RK4 vs ode45 Comparison: ^{99}Mo → ^{99m}Tc → ^{99}Tc');
subtitle('Transient Equilibrium Demonstration');
grid on;

%% Plot Number of Atoms (semilog)
figure;
% ode45 results 
semilogy(t_ode45, N_ode45(:,1), 'm-', 'LineWidth', 2, 'DisplayName', '^{99}Mo (ode45)'); hold on;
semilogy(t_ode45, N_ode45(:,2), 'c-', 'LineWidth', 2, 'DisplayName', '^{99m}Tc (ode45)');
semilogy(t_ode45, N_ode45(:,3), 'b-', 'LineWidth', 2, 'DisplayName', '^{99}Tc (ode45)');

% rk4 results 
semilogy(t_rk4, N_rk4(:,1), 'm--', 'LineWidth', 1.5, 'DisplayName', '^{99}Mo (RK4)');
semilogy(t_rk4, N_rk4(:,2), 'c--', 'LineWidth', 1.5, 'DisplayName', '^{99m}Tc (RK4)');
semilogy(t_rk4, N_rk4(:,3), 'b--', 'LineWidth', 1.5, 'DisplayName', '^{99}Tc (RK4)');
hold off;

xlabel('Time (hours)');
ylabel('Relative Number of Atoms');
legend('Location', 'best');
title('RK4 vs ode45 Comparison: Number of Atoms');
subtitle('Transient Equilibrium');
grid on;
