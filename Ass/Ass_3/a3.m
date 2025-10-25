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
N0 = [1/lambda_1, 0, 0];

% Time span (in days)
tSpan = [0, 1];

% Solve ODEs
[t, N] = ode45(decayChain, tSpan, N0);

% Compute activities (A = λ * N)
A_Ce = lambda_1 * N(:, 1);
A_Pr_m = lambda_2 * N(:, 2);
A_Pr = lambda_3 * N(:, 3);

figure;
plot(t, A_Ce, 'b', 'LineWidth', 1.5); hold on;
plot(t, A_Pr_m, 'r', 'LineWidth', 1.5);
plot(t, A_Pr, 'g', 'LineWidth', 1.5);
hold off;

legend('^{144}Ce', '^{144m}Pr', '^{144}Pr', 'Location', 'best');
xlabel('Time (days)');
ylabel('Activity (decays/day)');
title('Radioactive Decay Chain: ^{144}Ce → ^{144m}Pr → ^{144}Pr');
subtitle('Secular Equilibrium Condition');
grid on;
