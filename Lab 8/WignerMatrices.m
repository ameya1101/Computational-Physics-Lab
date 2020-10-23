M = 50000;
eigs = zeros(1, M);

subplot(2, 2, 1)
eigs1 = eigDistribution(5, M);
histogram(eigs1, 50);
title('$N=5$', 'interpreter', 'latex');
xlabel('Eigenvalues');
ylabel('Counts');

subplot(2, 2, 2)
eigs2 = eigDistribution(10, M);
histogram(eigs2, 50);
title('$N=10$', 'interpreter', 'latex');
xlabel('Eigenvalues');
ylabel('Counts');

subplot(2, 2, 3)
eigs3 = eigDistribution(20, M);
histogram(eigs3, 50);
title('$N=20$', 'interpreter', 'latex');
xlabel('Eigenvalues');
ylabel('Counts');

subplot(2, 2, 4)
eigs4 = eigDistribution(50, M);
histogram(eigs4, 50);
title('$N=50$', 'interpreter', 'latex');
xlabel('Eigenvalues');
ylabel('Counts');


function eigs = eigDistribution(N, M)
    tol = 1e-4;
    eigs = zeros(1, M);
    for i = 1:M
        A = rand(N);
        A_tran = transpose(A);
        S = (A + A_tran)/2;

        [eig, ~] = norm_power_method(S, tol);
        eigs(i) = eig;

    end
end

%hist(eigs);


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
            break;
        end
        
        eig0 = eig1;
        N = N + 1;
    end
    
    eig = eig1;
    eigvector = y;
   
end