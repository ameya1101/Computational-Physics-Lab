% Define the constants
m = 0.1;    % The mass of the object
h0 = 10;    % The height from which the object is dropped
b = 1.9e-4; % The damping constant

place = "Kuala-Lampur"; % Define the place at which t_fall has to be computed
g = find_g(place);

% Define the function
h = @(t) h0 - (m*g/b)*t + (m^2/b^2)*g*(1 - exp(- b*t/m));

% Ask user for two guesses to the root
guesses = input("Enter two initial guesses (within a []): ");
t0 = guesses(1);
t1 = guesses(2);

tol = 1e-4; % Minimum tolerance for termination

% Run the algorithm
[zero, flag] = SecantMethod(h, t0, t1, tol);


%Print the roots and plots
if flag
    fprintf("t_fall: %.4f\n", zero);
else
    fprintf("t_fall could not be found. The last computed guess was: %.4f\n", zero);
end

% Secant Method Algorithm
function [zero, flag] = SecantMethod(f, x0, x1, tol)
    N_max = 100; % Maximum number of iteration for termination
    N = 1;
    flag = 0;    % flag = 1 if a root is found, else 0
    
    while N < N_max
        
        f_x0 = f(x0);
        f_x1 = f(x1);
        
        x2 = x1 - f_x1*(x1 - x0)/(f_x1 - f_x0);
        
        error = abs(x2 - x1);
        
        if(error <= tol)
            zero = x2;
            fprintf("Number of iterations required for convergence: %d\n", N);
            fprintf("Error at termination: %f\n", error);
            flag = 1;
            break;
        end
        
        x0 = x1;
        x1 = x2;
       
    end
        
end

function g = find_g(place)
    if(nargin == 0 || place == "")
        g = 9.806; % Default to an average value of g if no argument is provided.
    else
        if(place == "Helsinki")
            g = 9.825; % Value of g in Helsinki, Finland
        elseif(place == "Toronto")
            g = 9.807; % Value of g in Toronto, Canada
        elseif(place == "Kuala-Lampur")
            g = 9.776; % Value of g in Kuala-Lampur, Malaysia
        end
    end
end