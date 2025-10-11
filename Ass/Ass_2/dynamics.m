Lambda = 5e-5;
beta = [0.00021,0.00142,0.00128,0.00257,0.00075,0.00027];
beta_total = sum(beta);
t_half = [56,23,6.2,2.3,0.61,0.23];
lambda_i = log(2)./t_half;
alpha = -4.2e-5;
Mfcf = 32.0;
tau = 4.5;
Rf = 0.1406;
P0 = 10.0;
Ti = 549.9997;
Tf0 = 551.4060;
rho_step = 0.2 * beta_total;   % convert 0.2 dollars to fraction
t_final = 1000;
output_dt = 50;
output_times = 0:output_dt:t_final;

Ci0 = (beta ./ (lambda_i * Lambda)) * P0;
dt_list = [10.0, 1.0, 0.5, 0.1];

for idt = 1:length(dt_list)
    dt = dt_list(idt);
    Nsteps = ceil(t_final/dt);
    times = (0:Nsteps)*dt;
    P = zeros(1, Nsteps+1);
    Tf = zeros(1, Nsteps+1);
    rho = zeros(1, Nsteps+1);
    Ci = zeros(6, Nsteps+1);

    P(1) = P0;
    Tf(1) = Tf0;
    Ci(:,1) = Ci0(:);
    rho(1) = rho_step + alpha*(Tf(1) - Ti);

    for n = 1:Nsteps
        A = (rho(n) - beta_total) / Lambda;
        SUMn = sum(lambda_i .* Ci(:,n)');

        % avoid overflow in exp
        argA = A*dt;
        if argA > 700
            argA = 700;
        elseif argA < -700
            argA = -700;
        end
        if abs(A) > 1e-14
            expA = exp(argA);
            P_particular = (SUMn / A) * (expA - 1.0);
            P(n+1) = P(n)*expA + P_particular;
        else
            P(n+1) = P(n) + SUMn * dt;
        end

        P_avg = 0.5*(P(n)+P(n+1));
        denom = (1.0 - rho(n));
        if abs(denom) < 1e-10
            denom = sign(denom)*1e-10;
            if denom == 0, denom = 1e-10; end
        end
        for i = 1:6
            expL = exp(-lambda_i(i)*dt);
            term = (1 - expL) / lambda_i(i);
            Ci(i,n+1) = Ci(i,n)*expL + (beta(i)/(Lambda*denom)) * P_avg * term;
        end

        expT = exp(-dt/tau);
        Tf(n+1) = Tf(n)*expT + ( (P(n)+P(n+1))/2 / Mfcf + Ti/tau ) * ( tau * (1 - expT) );
        rho(n+1) = rho(n) + alpha*(Tf(n+1) - Tf(n));
    end

    idx_out = arrayfun(@(t) find(abs(times - t) < 1e-8, 1), output_times);
    Time_out = times(idx_out)';
    Power_out = P(idx_out)';
    Tf_out = Tf(idx_out)';
    rho_out = rho(idx_out)' / beta_total;   % convert to dollars for output

    T = table(Time_out, Power_out, Tf_out, rho_out, ...
        'VariableNames', {'time_s','Power_W','FuelTemp_C','Reactivity_dollars'});

    filename = sprintf('sim_results_dt_%0.2fs.csv', dt);
    writetable(T, filename);
    fprintf('Saved results dt = %.2f s -> %s\n', dt, filename);
end

fprintf('\nSensitivity summary (final values at t = %d s):\n', t_final);
for idt = 1:length(dt_list)
    dt = dt_list(idt);
    folderfile = sprintf('sim_results_dt_%0.2fs.csv', dt);
    if exist(folderfile, 'file')
        TT = readtable(folderfile);
        finalRow = TT(end,:);
        fprintf(' dt=%.2fs -> P=%.6g W, Tf=%.6g °C, rho=%.6g $\n', ...
            dt, finalRow.Power_W, finalRow.FuelTemp_C, finalRow.Reactivity_dollars);
    end
end
