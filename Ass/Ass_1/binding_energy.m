function BE = binding_energy(Z, A, mass, isAtomicMass)
    mp = 1.007276466812;   % Proton mass in amu
    mn = 1.00866491595;    % Neutron mass in amu
    me = 0.000548579909;   % Electron mass in amu
    amu_to_MeV = 931.494;  % 1 amu = 931.494 MeV

    % number of neutrons
    N = A - Z;

    if isAtomicMass
        M_nucleus = mass - Z*me;  % Nuclear mass
    else
        M_nucleus = mass;   % atomic mass
    end

    mass_diff = (Z*mp + N*mn) - M_nucleus;

    % binding energy
    BE = mass_diff * amu_to_MeV;
end