%% Ask the user for data
    % User enters the function, the left interval limit and the right
    % interval limit. 
func_str = input("Enter the function: ", 's');
func = str2sym(func_str);
x_l = str2double(input("Enter x_L: ", 's'));
x_r = str2double(input("Enter x_R: ", 's'));

%% Run the algorithm
[zero, flag] = BisectionMethod(func, x_l, x_r);

%% Print the roots and plots
if flag
    fprintf("The root in the given interval: %f\n", zero);
else
    fprintf("The root could not be found in the given interval\n");
end

%% Bisection Method Algorithm
function [zero, flag] = BisectionMethod(f, x_l, x_r)
    epsilon_0 = 5e-10;
    delta_0 = 2e-10;    
    N_max = 100;
    N = 0;
    zero = 0;
    flag = 0;
    
    epsilons = [];
    iters = [];
    deltas = [];
    
    while N < N_max
        
        f_l = subs(f, x_l);
        f_r = subs(f, x_r);
        
        x_mid = (x_l + x_r)/2;
        f_mid = subs(f, x_mid);
        
        
        delta = abs(f_mid);
        epsilon = abs(x_l - x_r);
        
        if ((subs(f, x_mid) == 0) || (epsilon < epsilon_0 && delta < delta_0))
            zero = x_mid;
            flag = 1;
            fprintf("Number of iterations required to converge to a root: %d\n", N);
            plot_error(); % Plot the errors as a function of the number of iterations
            break;
            
        elseif f_mid * f_l < 0 
            x_r = x_mid;
            
        elseif f_mid * f_r < 0
            x_l = x_mid;
        end
        
        iters = [iters, N];
        epsilons = [epsilons, epsilon];
        deltas = [deltas, delta];
        
        N = N + 1;
    end
    
    
    function plot_error
        hold on
        plot(iters, epsilons)
        plot(iters, deltas)
        legend({'$\epsilon = |x_L - x_R|$', '$\delta = |f(x_{mid})|$'},'Interpreter','latex');
        ylabel("Errors", "interpreter", "latex");
        xlabel("Iterations $N$", "interpreter", "latex");
        hold off
    end
    
end




