%% Ask the user for data
    % User enters the function and an initial guess for the root.
func_str = input("Enter the function: ", 's');
func = str2sym("@(x)" + func_str);
x0 = str2double(input("Enter the initial guess: ", 's'));

%% Run the algorithm
[zero, flag] = NRMethod(func, x0);

%% Print the roots and plots
if flag
    fprintf("The root: %f\n", zero);
else
    fprintf("The root could not be found.\n");
end

%% Newton-Raphson Method Algorithm
function [zero, flag] = NRMethod(f, x0)
    epsilon = 5e-10;    
    N_max = 100;
    N = 1;
    zero = 0;
    flag = 0;
    
    errors = [];
    iters = [];
    
    x = x0;
    x_old = x0;
    
    f_diff = diff(f);
    
    while N < N_max
        
        f_x0 = subs(f, x);
        f_diff_x0 = subs(f_diff, x);
        
        x = x_old - f_x0/f_diff_x0;
        
        error = abs(x - x_old);
        
        if(error <= epsilon)
            zero = x;
            fprintf("Number of iterations required for convergence: %d\n", N);
            flag = 1;
            plot_error();
            break;
        end
        
        x_old = x;
       
        errors(N) = error;
        iters(N) = N;
        N = N + 1;
    end
    
    
    function plot_error
        hold on
        plot(iters, errors)
        legend('$\epsilon = |x_{i+1} - x_i|$','Interpreter','latex');
        ylabel("Error", "interpreter", "latex");
        xlabel("Iterations $N$", "interpreter", "latex");
        hold off
    end
    
end