A = [3 0 0; 0 -11 0; 0 0 4];
tol = 1e-4;

[eig, eigvector] = norm_power_method(A, tol);

fprintf("Matrix:  \n");
disp(A);
fprintf("Eigenvalue: %f\n", eig);
fprintf("Eigenvector: \n");
disp(eigvector);

function [eig, eigvector] = norm_power_method(A, tol)
    
    length = size(A);
    dim = length(1);          % Determine the dimensions of the matrix
    
    x = ones(dim, 1);        % Arbitary vector x0 to start computation
    y = zeros(dim, 1);
    eig0 = 42;

    N = 0;
    N_max = 1000;             % Maximum number of iterations

    while N < N_max
        y = A*x;
        eig1 = max(abs(y));
        x = y/eig1;
        
        if(abs(eig1 - eig0) < tol)
            fprintf("Iterations required for convergence: %d\n", N);
            break;
        end
        
        eig0 = eig1;
        N = N + 1;
    end
    
    eig = eig1;
    eigvector = y;
   
end