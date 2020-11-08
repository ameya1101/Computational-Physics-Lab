% Define initial and final simulation times. Step size is defined here,
% although it can be implemented inside the function itself.
t0 = 0;
tf = 70;
h = 0.01;

% Parameters for the Lorenz attractor.
r = 0.5; 
p = 10; 
b = 8/3;

% Lorenz differential equations and initial values.
f = @(x)[-p*x(1) + p*x(2); ...
          r*x(1) - x(2) - x(1)*x(3);...
         -b*x(3) + x(1)*x(2)];
inits = [1; 1; 1];

% Function call
[sol, t] = RK3(f, t0, tf, inits, h);

% Plot the required graphs as subplots for convenience.
subplot(1, 2, 1)
plot(t, sol(1, :));
title('Z vs. time', 'interpreter', 'latex');
xlabel('time', 'interpreter', 'latex');
ylabel('Z', 'interpreter', 'latex');

subplot(1, 2, 2)
plot(sol(1, :), sol(3, :))
title('Z vs. X', 'interpreter', 'latex');
xlabel('X', 'interpreter', 'latex');
ylabel('Z', 'interpreter', 'latex');
box on


% Define the function for computing the RK3 algorithm

function [sol, t] = RK3(f, t0, tf, inits, h)
% RK3 Find the solution to a system of ordinary differential equations
% using the Runge-Kutta 3 method.
%
% 	[sol, t] = RK(f, t0, tf, inits, stepsize) takes five mandatory parameter
% 	values.
%       f -> the array of functions 
%       t0 -> the start time of the simulation
%       tf -> the end time of the simulation
%       inits -> the initial value of the functions at t0
%       h -> the stepsize for the algorithm
%
%   This function implements Heun's third-order method given by the
%   Butcher tableau
%   
%   0   |  0   0   0
%   1/3 | 1/3  0   0 
%   2/3 |  0  2/3  0
%   -------------------
%       | 1/4  0  3/4
%   

    
    % Determine the number of functions
    sz = size(f(inits));
    N = sz(1);                              % Number of ODEs in the system
    
    % Determine the number of steps.
    steps = ceil((tf - t0)/h);
    
    % Create solution arrays and initialise them.
    t = linspace(t0, tf, steps);
    sol = zeros(N, steps);
    sol(:, 1) = inits;
    
    % Perform the iterative steps.
    for step = 1 : steps - 1
        k1 = f(sol(:, step));
        k2 = f(sol(:, step) + 1/3 * k1 * h);
        k3 = f(sol(:, step) + 2/3 * k2 * h);
        sol(:, step+1) = sol(:, step) + h * (1/4 * k1 + 3/4 * k3);
    end
    
end

