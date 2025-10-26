function [t, y] = rk4_solver(odefun, tspan, y0, h)
% rk4_solver Solves a system of ODEs using the classic 4th-order Runge-Kutta method.
%
%   Inputs:
%       odefun - Function handle for the ODE system, f(t, y).
%       tspan  - Time interval [t0, tf].
%       y0     - Initial conditions (column vector).
%       h      - Step size (scalar).
%
%   Outputs:
%       t      - Vector of time points.
%       y      - Solution matrix, where each row corresponds to a time
%                point in t and each column is a state variable.

    % Ensure y0 is a column vector for calculations
    if isrow(y0)
        y0 = y0';
    end

    % Create the time vector
    t0 = tspan(1);
    tf = tspan(2);
    t = (t0:h:tf)';
    
    % Get number of steps and number of equations
    num_steps = length(t);
    num_eqns = length(y0);
    
    % Initialize the solution matrix
    y = zeros(num_steps, num_eqns);
    
    % Set the initial condition
    y(1, :) = y0'; % Store as a row
    
    % Perform the RK4 iteration
    for i = 1:(num_steps - 1)
        % Get current time and state (as a column vector)
        ti = t(i);
        yi = y(i, :)';
        
        % Calculate the four RK4 "slopes"
        k1 = odefun(ti, yi);
        k2 = odefun(ti + 0.5*h, yi + 0.5*h*k1);
        k3 = odefun(ti + 0.5*h, yi + 0.5*h*k2);
        k4 = odefun(ti + h, yi + h*k3);
        
        % Calculate the next state
        y_next = yi + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
        
        % Store the next state as a row
        y(i+1, :) = y_next';
    end
end