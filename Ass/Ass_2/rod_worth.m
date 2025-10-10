% given data
%       x1   x2   t
data = [0     89    72.00;
        89   139   124.00;
        139  192    75.00;
        192  229   122.30;
        229  288    49.00;
        288  339    56.00;
        339  380    73.20;
        380  416   113.50;
        416  443   249.80;
        443  481   103.60;
        481  519   126.40;
        519  582    81.50;
        582  656    88.00;
        656  729   143.10;
        729  999   123.50;
        ];

beta_eff = 0.007;  % Effective delayed neutron fraction
ell = 3.2e-5;  % Prompt neutron lifetime (s)

beta = [0.000231, 0.00153, 0.001372, 0.00276, 0.000805, 0.000294];
lambda = [0.0124, 0.0305, 0.1115, 0.301, 1.138, 3.01];

n = size(data, 1);
rho_dollar = zeros(n, 1);
diff_worth = zeros(n, 1);
int_worth = zeros(n, 1);
pcm_conv = beta_eff * 1e5;

for i = 1:n
    pos1 = data(i, 1);
    pos2 = data(i, 2);
    t = data(i, 3);

    % Reactor period T
    T = t / log(4); % P/P0 = 800/200 = 4

    % Inhour equation
    sum_term = 0;
    for j = 1:6
        sum_term = sum_term + beta(j) / (1 + lambda(j)*T);
    end

    rho_dollar(i) = (1 / beta_eff) * ((ell/T + sum_term) / (1 + ell/T));

    % Differential worth per step
    dx = pos2 - pos1;
    diff_worth(i) = rho_dollar(i) / dx;

    % Integral worth 
    if i == 1
        int_worth(i) = rho_dollar(i);
    else
        int_worth(i) = int_worth(i-1) + rho_dollar(i);
    end
end

diff_worth_pcm = diff_worth * pcm_conv;
int_worth_pcm = int_worth * pcm_conv;

disp(' Stroke  dx      t(s)      rho($)    Diff($/step)   Int($)');
for i = 1:n
    fprintf('%3d    %4.0f    %7.2f    %8.5f    %8.5f    %8.5f\n', ...
        i, data(i, 2)-data(i, 1), data(i, 3), rho_dollar(i), diff_worth(i), int_worth(i));
end

fprintf('\nTotal Integral Worth = %.5f $ = %.2f pcm\n', int_worth(end), int_worth(end)*pcm_conv);
