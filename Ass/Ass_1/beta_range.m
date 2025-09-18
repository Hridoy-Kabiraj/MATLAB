% Approximate density of mediums
ro_air = 0.001225; %g/cm3
ro_polyethylene = 0.93; %g/cm3
ro_aluminium = 2.70; %g/cm3
ro_iron = 7.87; %g/cm3
ro_lead = 11.34; %g/cm3

max_energy = input("Enetr the maximum enery of beta particle: ");
energy = 0.01 : 0.00001 : max_energy;
N = length(energy);

air_range = zeros(1, N);
polyethylene_range = zeros(1, N);
aluminium_range = zeros(1, N);
iron_range = zeros(1, N);
lead_range = zeros(1, N);


for i = 1:N
    E = energy(i);
    if E <= 2.5
        range = (412*E^(1.265 - 0.0954*log(E)))*10^(-3);  %g/cm2
    else
        range = (530*E -106)*10^(-3); %g/cm2
    end
    air_range(i) = range / ro_air;  %cm
    polyethylene_range(i) = range / ro_polyethylene; %cm
    aluminium_range(i) = range / ro_aluminium; %cm
    iron_range(i) = range / ro_iron;  % cm
    lead_range(i) = range / ro_lead; %cm
end

figure;
semilogy(energy, air_range,'r', 'LineWidth', 1.5);
hold on;
semilogy(energy, polyethylene_range, 'b ', 'LineWidth', 1.5);
semilogy(energy, aluminium_range, 'g', 'LineWidth', 1.5);
semilogy(energy, iron_range, 'k', 'LineWidth', 1.5);
semilogy(energy, lead_range, 'm ', 'LineWidth', 1.5);
xlabel('Beta Energy (MeV)');
ylabel('Range(cm)');
title('Beta Particle Range in Different Materials');
legend('Air', 'Ployethylene', 'Aluminium', 'Iron', 'Lead');
grid on;

    