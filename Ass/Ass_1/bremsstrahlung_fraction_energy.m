% effective atomic number for materials
z_air = 7.6; %approximate
z_polyethylene = 5.5; %approximate
z_aluminium = 13;
z_iron = 26;
z_lead = 82;

max_energy = input("Enter the maximum energy of beta particle in MeV: "); %MeV
energy = 0.01 : 0.00001 : max_energy;
N = length(energy);

% initializing Bremsstrahlung fraction-energy with zeros
brem_air = zeros(1, N);
brem_polyethylene = zeros(1, N);
brem_aluminium = zeros(1, N);
brem_iron = zeros(1, N);
brem_lead = zeros(1, N);

for i = 1:N
    E = energy(i);
    brem_air(i) = (E*z_air) / 800;
    brem_polyethylene(i) = (E*z_polyethylene) / 800;
    brem_aluminium(i) = (E*z_aluminium) /800;
    brem_iron(i) = (E*z_iron) / 800;
    brem_lead(i) = (E*z_lead) / 800;
end

figure;
plot(energy, brem_air, 'r', 'LineWidth', 1.5);
hold on;
plot(energy, brem_polyethylene, 'b ', 'LineWidth', 1.5);
plot(energy, brem_aluminium, 'g', 'LineWidth', 1.5);
plot(energy, brem_iron, 'k', 'LineWidth', 1.5);
plot(energy, brem_lead, 'm ', 'LineWidth', 1.5);
xlabel('Beta Particle Enenrgy(MeV)');
ylabel('Bremsstrahlung fraction-energy');
title('Bremsstrahlung fraction-energy curve for beta particles');
legend('Air', 'Polyethylene', 'Aluminium', 'Iron', 'Lead', 'Location','best');
grid on;