% Ask user for a function and two guesses for the root
func = str2sym("@(x)" + input("Enter the function: ", 's'));
guesses = input("Enter two initial guesses (within a []): ");
x0 = guesses(1);
x1 = guesses(2);

% Run the algorithm
[zero, flag] = SecantMethod(func, x0, x1);

% Print the roots and plots
if flag
    fprintf("The root: %f\n", zero);
else
    fprintf("The root could not be found.\n");
end

% Newton-Raphson Method Algorithm
function [zero, flag] = SecantMethod(f, x0, x1)
    epsilon = 5e-10;    
    N_max = 100;
    N = 1;
    zero = 0;
    flag = 0;
    
    errors = [];
    iters = [];
    
    x0 = x0;
    x1 = x1;
    x2 = 0;
    
    
    while N < N_max
        
        f_x0 = subs(f, x0);
        f_x1 = subs(f, x1);
        f_x2 = subs(f, x2);
        
        x2 = x1 - f_x1*(x1 - x0)/(f_x1 - f_x0);
        
        error = abs(x2 - x1);
        
        if(error <= epsilon)
            zero = x2;
            fprintf("Number of iterations required for convergence: %d\n", N);
            flag = 1;
            plot_error();
            break;
        end
        
        x0 = x1;
        x1 = x2;
       
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