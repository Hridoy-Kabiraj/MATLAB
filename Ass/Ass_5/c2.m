%% Constants and Definitions
% Homogenized Number Densities ( from 2a ) , in atoms / cm ^3
N_C12 = 8.40e22 ; % Carbon -12 ( Graphite )
N_U235 = 3.82e18 ; % Uranium -235
N_U238 = 5.26e20 ; % Uranium -238
barn_to_cm2 = 1e-24; % 1 barn = 10^ -24 cm ^2

%% Load Microscopic Cross Section Data
c12_elastic_data = readmatrix ( 'C12_elastic.csv' );
u235_capture_data = readmatrix ( 'U235_capture.csv' );
u235_fission_data = readmatrix ( 'U235_fission.csv' );
u238_capture_data = readmatrix ( 'U238_capture.csv' );
u238_fission_data = readmatrix ( 'U238_fission.csv' );
u238_elastic_data = readmatrix ( 'U238_elastic.csv' );

%% Calculate Macroscopic Cross Sections ( Sigma )
E_c12_el = c12_elastic_data (: , 1) ;
sig_c12_el = c12_elastic_data (: , 2) ;
SIG_c12_el = N_C12 * sig_c12_el * barn_to_cm2 ;
E_u235_cap = u235_capture_data (: , 1) ;
sig_u235_cap = u235_capture_data (: , 2) ;
SIG_u235_cap = N_U235 * sig_u235_cap * barn_to_cm2 ;
E_u235_fis = u235_fission_data (: , 1) ;
sig_u235_fis = u235_fission_data (: , 2) ;
SIG_u235_fis = N_U235 * sig_u235_fis * barn_to_cm2 ;
E_u238_cap = u238_capture_data (: , 1) ;
sig_u238_cap = u238_capture_data (: , 2) ;
SIG_u238_cap = N_U238 * sig_u238_cap * barn_to_cm2 ;
E_u238_fis = u238_fission_data (: , 1) ;
sig_u238_fis = u238_fission_data (: , 2) ;
SIG_u238_fis = N_U238 * sig_u238_fis * barn_to_cm2 ;
E_u238_el = u238_elastic_data (: , 1) ;
sig_u238_el = u238_elastic_data (: , 2) ;
SIG_u238_el = N_U238 * sig_u238_el * barn_to_cm2 ;

%% Calculate Neutron Energy Spectra ( for 2c )
E_spec = logspace ( -3 , 7.3 , 10000) ; % From 0.001 eV to ~20 MeV
% 1. Maxwell - Boltzmann Spectrum at T =300 K
k = 8.617333e-5; % Boltzmann constant in eV / K
T = 300; % Temperature in Kelvin
kT = k * T ; % Thermal energy in eV
phi_MB = E_spec .* exp ( - E_spec / kT ) ; % Unscaled Maxwellian flux
% 2. Watt Fission Spectrum for U -235
E_MeV = E_spec * 1e-6; % Convert energy grid to MeV for Watt formula
a = 0.983; % Constant a in MeV ^ -1
b = 2.29; % Constant b in MeV ^ -1
chi_Watt = exp ( -a * E_MeV ) .* sinh ( sqrt ( b * E_MeV ) ) ; % Unscaled Watt spectrum
% Scale the spectra to be visually informative on the plot
max_Sigma = max ( SIG_u235_fis ) ;
phi_MB_scaled = phi_MB * ( max_Sigma / max ( phi_MB ) ) * 0.5;
chi_Watt_scaled = chi_Watt * ( max_Sigma / max ( chi_Watt ) ) * 0.8;

%% Plotting the Results
figure ( 'Name' , 'Macroscopic Cross Sections and Spectra for CP-1' );
hold on ;
loglog ( E_c12_el , SIG_c12_el , 'k-' , 'LineWidth' , 2 , 'DisplayName' , 'C-12 Elastic' );
loglog ( E_u235_cap , SIG_u235_cap , 'b--' , 'LineWidth' , 1.5 , 'DisplayName' , 'U-235 Capture' );
loglog ( E_u235_fis , SIG_u235_fis , 'b-' , 'LineWidth' , 2 , 'DisplayName' , 'U-235 Fission' );
loglog ( E_u238_cap , SIG_u238_cap , 'r--' , 'LineWidth' , 1.5 , 'DisplayName' , 'U-238 Capture' );
loglog ( E_u238_fis , SIG_u238_fis , 'r-' , 'LineWidth' , 2 , 'DisplayName' , 'U-238 Fission' );
loglog ( E_u238_el , SIG_u238_el , 'r:' , 'LineWidth' , 1.5 , 'DisplayName' , 'U-238 Elastic' );
loglog ( E_spec , phi_MB_scaled , 'g-' , 'LineWidth' , 2.5 , 'DisplayName' , 'Maxwellian Spectrum (scaled)' );
loglog ( E_spec , chi_Watt_scaled , 'm-' , 'LineWidth' , 2.5 , 'DisplayName' , 'Watt Spectrum (scaled)' );
hold off ;
title ( 'Macroscopic Cross Sections & Spectra for Homogenized CP-1 Mixture' );
xlabel ( 'Neutron Energy (eV)' );
ylabel ( 'Macroscopic Cross Section \Sigma (cm^{-1}) / Scaled Flux' );
legend ( 'show' , 'Location' , 'southwest' );
grid on ;
set ( gca , 'XScale' , 'log' , 'YScale' , 'log' );
xlim ([1e-3 2e7 ]) ;
disp ( 'Plot generated successfully.' );
