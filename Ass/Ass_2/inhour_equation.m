Lambda = 1e-5;  % Prompt neutron generation time
beta = [0.000266, 0.001491, 0.001316, 0.002849, 0.000896, 0.000182]; % Delayed neutron fraction
lambda = [0.0127, 0.0317, 0.115, 0.311, 1.4, 3.87]; % Decay constant

% The reactor root
s = logspace(-3, 2, 500); % vary over 10^-3 to 10^2

% Inhour equation
rho = s.*Lambda;
for i = 1:6
    rho = rho +  s.*(beta(i)./(s + lambda(i)));
end

rho_pcm = rho * 1e5;

% Plot
figure
semilogx(1./s, rho_pcm, 'Linewidth', 2)
xlabel('Ractor Period T = 1/s (s)')
ylabel('Reactivity (pcm)')
title('Inhour Equation')
grid on