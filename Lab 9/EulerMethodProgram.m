N = 2;

g = 9.81;
L = 1;


funcs = cell(N, 1);
funcs{1} = @(t, x) x;
funcs{2} = @(t, x) - (g/L)*x;

inits = [0, 1];
h = 0.01;
t_max = 2;

[v, t] = EulerMethod(funcs, N, inits, t_max, h);

hold on
title("$f(t) = e^{-t}$", 'interpreter', 'latex');
xlabel('time');
ylabel('$f(t)$', 'interpreter', 'latex');
plot(t, y(2, :));
hold off


function [y, t] = EulerMethod(funcs, N, inits, t_max, step_size)
    y = zeros(N, t_max/step_size);
        for i = 1:N
            t = 0;
            ctr = 0;
            while (t < t_max - step_size)
               y(i, ctr+1) = inits(i) + step_size .* funcs{i}(inits(i), ctr);
               inits(i) = y(i, ctr+1);
               
               ctr = ctr + 1;
               t = t + step_size;
            end
        end
      t = linspace(0, t_max, t_max/step_size);
end