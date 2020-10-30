g = 9.81;
L = 1;

y0 = [0; 1];
t_max = 100;
stepsize = 0.1;

funcs = @(x)[-(g/L)*x(2); x(1)];

[y, t] = RK2(funcs, t_max, y0, stepsize);
hold on
subplot(1, 2, 1)
plot(t, y(2, :));
title('Displancement vs Time', 'interpreter', 'latex');
xlabel('Time', 'interpreter', 'latex');
ylabel('Displacement', 'interpreter', 'latex');
subplot(1, 2, 2)
plot(y(1, :), y(2, :));
title('Phase Space', 'interpreter', 'latex');
xlabel('Velocity', 'interpreter', 'latex');
ylabel('Displacement', 'interpreter', 'latex');
hold off

function [y, t] = RK2(funcs, t_max, y0, stepsize)

    sz = size(funcs(y0));
    N = sz(1);
    
    steps = ceil(t_max/stepsize);
    t = zeros(steps, 1);
    y = zeros(N, steps);
    y(:, 1) = y0;
    
    for step = 1 : steps - 1
        k1 = funcs(y0);
        k2 = funcs(y0 + k1 * stepsize/2);
        y(:, step) = y0 + stepsize * k2;
        y0 = y(:, step);
        t(step+1) = t(step) + stepsize;
    end
    
end

