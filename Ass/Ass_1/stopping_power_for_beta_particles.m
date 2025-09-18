% Defining the parameters
r_e = 2.8179e-13; %cm  classical radius of electron
m_e = 9.109e-28; %g  mass of electron
c = 2.998e10; %cm/s velocity of light
z = 1;  % atomic number of the projectile (beta)
Z = 45.798; % target effective atomic number (NaI)
I = 452e-6; %MeV ionization potential
me_c2 = 0.511; %Mev rest energy of electron

% Calculating N
rho = 3.67; %g/cm3 for NaI
N_A = 6.023e23;
A = 149.89; % kg/mol molar mass
N = rho * N_A / A;

erg_to_MeV = 6.24151e5; % 1 erg = 6.24151×10⁵ MeV

function v = velocity(E, c, me_c2)
    KE = E; %Kinetic energy of beta particle

    gamma = 1.0 + KE / me_c2;
    v = c * sqrt(1 - (1 / gamma^2));
end

max_energy = input("Enter the maximum energy of Beta particle in MeV: ");
energy = 1.0 : 0.1 : max_energy;
len = length(energy);
S_e = zeros(1, len); % The stopping powe r for electrons is due to the electronic ionization
S_r = zeros(1, len); % The stopping power for electrons is due to the electronic radiative loss

for i = 1:len
    E = energy(i);
    v = velocity(E, c, me_c2);
    beta = v / c;
    beta2 = beta^2;
    beta1 = 1 - beta2;

    cons = (2*pi*r_e^2*me_c2*N*z^2*Z) / beta2;
    mv2_MeV = (m_e * v^2) * erg_to_MeV;
    log_term = log((mv2_MeV * E) / (2 * I^2 * beta1));
    correction_term = -log(2) * (beta2 - 1 + 2*sqrt(1-beta2));
    additional_term = (1/8) * (1 - sqrt(1-beta2));
    S_e(i) = cons * (log_term + correction_term + beta1 + additional_term);

    rad_const = ((z+1) * z * rho * E) / (1847 * me_c2^2);
    rad_log_term = 4 * log(2*E/me_c2) - 4/3;
    S_r(i) = rad_const * rad_log_term;
end

figure;
plot(energy, S_e, 'r', 'LineWidth', 1.5);
hold on;
plot(energy, S_r, 'b ', 'LineWidth', 1.5);
xlabel('Energy (MeV)');
ylabel('Stopping Power (MeV/cm)');
legend('Electronic','Radiative');
grid on;
title('Stopping Power of Beta Particles in NaI');