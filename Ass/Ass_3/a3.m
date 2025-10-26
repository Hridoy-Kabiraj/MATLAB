clc; clearvars; close all;

% Half-lives (in days)
T_Ce = 284.893;                  % Ce-144
T_Pr_m = 7.2 / (24 * 60);        % Pr-144m (7.2 hours)
T_Pr = 17.29 / (24 * 60);        % Pr-144 (17.29 minutes)

% Decay constants (per day)
lambda_1 = log(2) / T_Ce;
lambda_2 = log(2) / T_Pr_m;
lambda_3 = log(2) / T_Pr;

% ODE system
decayChain = @(t, N) [
    -lambda_1 * N(1);                                 % Ce-144
     0.0115 * lambda_1 * N(1) - lambda_2 * N(2);      % Pr-144m
     0.9885 * lambda_1 * N(1) - lambda_3 * N(3)       % Pr-144
];

% Initial condition (atoms)
N0 = [1/lambda_1; 0; 0]; % Use column vector for consistency

% Time span (in days)
tSpan = [0, 1];

% Solve with ode45
[t_ode45, N_ode45] = ode45(decayChain, tSpan, N0);
% Compute activities
A_Ce_ode45 = lambda_1 * N_ode45(:, 1);
A_Pr_m_ode45 = lambda_2 * N_ode45(:, 2);
A_Pr_ode45 = lambda_3 * N_ode45(:, 3);

% Solve with rk4_solve
h = 0.001; % Define a step size for RK4
[t_rk4, N_rk4] = rk4_solver(decayChain, tSpan, N0, h);
% Compute activities
A_Ce_rk4 = lambda_1 * N_rk4(:, 1);
A_Pr_m_rk4 = lambda_2 * N_rk4(:, 2);
A_Pr_rk4 = lambda_3 * N_rk4(:, 3);


figure;
% Plot ode45 results 
plot(t_ode45, A_Ce_ode45, 'b-', 'LineWidth', 2, 'DisplayName', '^{144}Ce (ode45)'); hold on;
plot(t_ode45, A_Pr_m_ode45, 'r-', 'LineWidth', 2, 'DisplayName', '^{144m}Pr (ode45)');
plot(t_ode45, A_Pr_ode45, 'g-', 'LineWidth', 2, 'DisplayName', '^{144}Pr (ode45)');

% Plot rk4 results (
plot(t_rk4, A_Ce_rk4, 'b--', 'LineWidth', 1.5, 'DisplayName', '^{144}Ce (RK4)');
plot(t_rk4, A_Pr_m_rk4, 'r--', 'LineWidth', 1.5, 'DisplayName', '^{144m}Pr (RK4)');
plot(t_rk4, A_Pr_rk4, 'g--', 'LineWidth', 1.5, 'DisplayName', '^{144}Pr (RK4)');
hold off;

legend('Location', 'best');
xlabel('Time (days)');
ylabel('Activity (decays/day)');
title('RK4 vs ode45 Comparison: ^{144}Ce → ^{144m}Pr → ^{144}Pr');
subtitle('Secular Equilibrium Condition');
grid on;