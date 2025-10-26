% Number Densities ( atoms / cm ^3)
N_C12 = 8.40e22 ;
N_U235 = 3.82e18 ;
N_U238 = 5.26e20 ;
barn_to_cm2 = 1e-24;

% Load Data
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

%% Plotting
figure ( 'Name' , 'Macroscopic Cross Sections for CP-1' );
hold on ;
loglog ( E_c12_el , SIG_c12_el , 'k-' , 'LineWidth' , 2 , 'DisplayName' , 'C-12 Elastic' );
loglog ( E_u235_cap , SIG_u235_cap , 'b--' , 'LineWidth' , 1.5 , 'DisplayName' , 'U-235 Capture' );
loglog ( E_u235_fis , SIG_u235_fis , 'b-' , 'LineWidth' , 2 , 'DisplayName' , 'U-235 Fission' );
loglog ( E_u238_cap , SIG_u238_cap , 'r--' , 'LineWidth' , 1.5 , 'DisplayName' , 'U-238 Capture' );
loglog ( E_u238_fis , SIG_u238_fis , 'r-' , 'LineWidth' , 2 , 'DisplayName' , 'U-238 Fission' );
loglog ( E_u238_el , SIG_u238_el , 'r:' , 'LineWidth' , 1.5 , 'DisplayName' , 'U-238 Elastic' );
hold off ;
title ( 'Macroscopic Cross Sections for Homogenized CP-1 Mixture' );
xlabel ( 'Neutron Energy (eV)' );
ylabel ( 'Macroscopic Cross Section \Sigma (cm^{-1})' );
legend ( 'show' , 'Location' , 'southwest' );
grid on ;
set ( gca , 'XScale' , 'log' , 'YScale' , 'log' );
xlim ([1e-3 2e7 ]) ;
