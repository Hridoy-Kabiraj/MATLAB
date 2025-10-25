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
N0 = [1/lambda_Mo99, 0, 0];   % Initial number of nuclei

%% Time span
tSpan = linspace(0, 100, 1000);  % hours

%% Solve the ODE system
[t, N] = ode45(decayChain, tSpan, N0);

%% Plot Activity vs Time (linear)
figure;
plot(t, lambda_Mo99*N(:,1), 'm', ...
     t, lambda_Tc99m*N(:,2), 'c', ...
     t, lambda_Tc99*N(:,3), 'b', 'LineWidth', 1.5);

xlabel('Time (hours)');
ylabel('Activity (decays/hour)');
legend('^{99}Mo (Parent)', '^{99m}Tc (Intermediate)', '^{99}Tc (Daughter)');
title('Radioactive Decay Chain: ^{99}Mo → ^{99m}Tc → ^{99}Tc');
subtitle('Transient Equilibrium Demonstration');
grid on;

%% Plot Number of Atoms (semilog)
figure;
semilogy(t, N(:,1), 'm', 'DisplayName', '^{99}Mo'); hold on;
semilogy(t, N(:,2), 'c', 'DisplayName', '^{99m}Tc');
semilogy(t, N(:,3), 'b', 'DisplayName', '^{99}Tc');

xlabel('Time (hours)');
ylabel('Relative Number of Atoms');
legend('show');
title('Radioactive Decay Chain of ^{99}Mo to ^{99}Tc');
subtitle('Transient Equilibrium');
grid on;
