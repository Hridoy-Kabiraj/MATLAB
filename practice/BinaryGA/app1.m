clc;
clear;
close all;


%% Problem Definition

problem.CostFunction = @(x) MinOne(x);
problem.nVar = 100;


%% GA Parameters

params.nPop = 50;
params.MaxIt = 100;

params.pC = 1;    % percentage of children(offsprings)
params.mu = 0.02;
params.beta = 1;




%% Run GA

out = RunGA(problem, params);


%% Results

figure;
plot(out.bestcost, LineWidth=2);
xlabel('Iterations');
ylabel('Best Cost');
grid on;