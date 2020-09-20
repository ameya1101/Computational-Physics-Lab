% Coefficient and constant matrices
A = [95 -25 -50 0 0 0;
    -25 76 -1 -30 0 0;
    -50 -1 126 -55 0 0;
    0 -30 -55 200 -25 -50;
    0 0 0 -25 56 -1;
    0 0 0 -50 -1 106];

B = [10; 0; 0; 0; 0; 0];

%tolerance
tol = 1e-6;

sol = jacobi(A, B, tol);

fprintf("The solution to the linear system is: \n");
disp(sol);

%Function to find the solution
function sol = jacobi(A, B, tol)
    S = size(B);
    dim = S(1);
    
    N_max = 250;
    
    D = diag(diag(A));
    T = -D\(A - D);
    C = D\B;
    
    flag = checkDiagDominance(T);
    
    if(flag)
        fprintf("Matrix is not diagonally dominant. Please try again\n");
    else
        
        X0 = rand(dim, 1);
        
        for i=1:N_max
            
            X = C + T*X0;
            
            if (max(abs((X - X0))) < tol)
                fprintf("Number of iterations required for convergence: %d\n", i);
                break;
            end
 
            X0 = X;
        end
        
        sol = X;
    end
end


function flag = checkDiagDominance(A)
    flag = eigs(A, 1) >= 1;
end
