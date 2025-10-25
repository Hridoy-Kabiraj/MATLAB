
clc; clearvars; close all;

% Grid
Nx = 51; Ny = 51;
dx = 1/(Nx-1);
dy = dx;    % square grid

% Initialize
phi = zeros(Nx,Ny);
phi(1,:) = 1;    % top boundary y=1
% other boundaries already zero

% Convergence
tol = 1e-6;
err = Inf;
maxIter = 1e6;
iter = 0;

% Gauss-Seidel in-place
while err > tol && iter < maxIter
    err = 0;
    iter = iter + 1;
    % update interior points (i = rows, j = cols)
    for i = 2:Nx-1
        for j = 2:Ny-1
            old = phi(i,j);
            phi(i,j) = 0.25*(phi(i+1,j) + phi(i-1,j) + phi(i,j+1) + phi(i,j-1));
            err = max(err, abs(phi(i,j)-old));
        end
    end
end

fprintf('Steady solved in %d iterations, final max-change = %.3e\n', iter, err);

% Plot
x = linspace(0,1,Nx);
y = 1 - linspace(0,1,Ny);  % flip for display
[X,Y] = meshgrid(x,y);
figure;
contourf(X,Y,phi,20,'LineColor','none');
colorbar; xlabel('x'); ylabel('y');
title('Steady-state \xi(x,y)');
