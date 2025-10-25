clc; clearvars; close all;

% grid
Nx = 51; Ny = 51;
dx = 1/(Nx-1);

% time step
dt = 0.25 * dx^2;
alpha = dt/dx^2;

% initial and BCs
phi = zeros(Nx,Ny);
phi(1,:) = 1;      % top row = 1
phi_new = phi;

% convergence
tol = 1e-6;
change = Inf;
step = 0;
maxSteps = 5e6;

% sample error every 100 steps for plotting
sample_err = [];
sample_time = [];

while change > tol && step < maxSteps
    step = step + 1;
    % explicit update (use current phi)
    phi_new(2:end-1,2:end-1) = phi(2:end-1,2:end-1) + ...
        alpha*( phi(3:end,2:end-1) + phi(1:end-2,2:end-1) + ...
                phi(2:end-1,3:end) + phi(2:end-1,1:end-2) - ...
                4*phi(2:end-1,2:end-1) );
    change = max(abs(phi_new(:)-phi(:)));
    phi = phi_new;
    if mod(step,100)==0
        sample_err(end+1)=change;
        sample_time(end+1)=step*dt;
    end
end

fprintf('Transient reached tol after %d steps (time = %.5f s), final change = %.2e\n', ...
    step, step*dt, change);

% plot final field
x = linspace(0,1,Nx);
y = 1 - linspace(0,1,Ny);
[X,Y] = meshgrid(x,y);
figure; contourf(X,Y,phi,20,'LineColor','none'); colorbar;
xlabel('x'); ylabel('y'); title('Transient final \xi(x,y)');

% plot error decay (sampled)
figure; plot(sample_time, sample_err, '-o'); grid on;
xlabel('Time (s)'); ylabel('max change'); title('Transient convergence');
