h = 1e-3;
L = 1;
T = 1;

N = L/h;

temp = massDensity(L, h, N-1);
mu_bottom = [temp(2:end), Inf];
mu_top = massDensity(L, h, N-1);
A = diag(T./(h^2*mu_top),1) + diag(-2*T./(h^2*massDensity(L, h, N))) + diag(T./(h^2*mu_bottom), -1);

tol = 1e-6;

[eig, y] = inv_power_method(A, tol, 250.3424);

fprintf("Fundamental eigenfrequency: %f\n", sqrt(-eig));
  
hold on

plot(y);
title("Normal Mode");
xlabel("$x$", 'interpreter', 'latex');
ylabel("$y(x)$", 'interpreter', 'latex');

hold off

function [eig, eigvector] = inv_power_method(A, tol, guess)
    
    length = size(A);
    dim = length(1);          % Determine the dimensions of the matrix
    A_inv = inv(A - guess*eye(dim));
    
    x0 = ones(dim, 1);      % Arbitary vector x0 to start computation
    y_t = rand(dim, 1).';   % Arbitraty vector to extract the eigenvalue
    eig0 = 42;

    N = 0;
    N_max = 1000;           % Maximum number of iterations

    while N < N_max
        x_n = A_inv*x0;
        x = x0;
        x0 = x_n;

        eig1 = (y_t*x_n)/(y_t*x);

        if(abs(eig1 - eig0) < tol)
            fprintf("Iterations required for convergence: %d\n", N);
            break;
        end

        eig0 = eig1;

        N = N + 1;
    end
    
    eig = 1/eig1;
    eigvector = x_n;
   
end

function mu = massDensity(L, h, N)
    
    mu0 = 0.954e-3;
    delta = 0.5e-3;
    mu = zeros(1, N);
    
    for i=1:N
        mu(i) = mu0 + (i*h - L/2)*delta;
    end
end