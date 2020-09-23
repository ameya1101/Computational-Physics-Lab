A = [3 0 0; 0 -11 0; 0 0 4];
tol = 1e-4;

[eig, eigvector] = inv_power_method(A, tol);

fprintf("Matrix:  \n");
disp(A);
fprintf("Eigenvalue: %f\n", eig);
fprintf("Eigenvector: \n");
disp(eigvector);

function [eig, eigvector] = inv_power_method(A, tol)
    
    length = size(A);
    dim = length(1);          % Determine the dimensions of the matrix
    A_inv = inv(A);
    
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

