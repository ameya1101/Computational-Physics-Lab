% Constants
m = 9.101E-31; %mass of electron
hbar = 1.054571817E-34;
e = 1.602E-19; %charge of electron

%Variables
V_0 = 100; %Potential
a = 3E-10; % Size of the potential well

%the coeff for a*alpha and a*beta
k = ( sqrt( 2*m*e )*a ) / hbar;

%E and V_0 are in terms of eV
a_alpha = @(E) k*sqrt(E);	
a_beta = @(E) k*sqrt(V_0 - E);

%bounds of E
bounds = @(E) 0 <= E && E <= V_0;

%functions in terms of a*alpha and a*beta
f_even = @(E) a_alpha(E)*tan( a_alpha(E) ) - a_beta(E);
f_odd = @(E) a_alpha(E)*cot( a_alpha(E) ) - a_beta(E);

%derivatives
df_even = @(E) ( (k^2)*0.5 ) * ( sec( a_alpha(E) )^2 + tan( a_alpha(E) )/a_alpha(E) + 1/a_beta(E) );
df_odd = @(E) ( (k^2)*0.5 ) * ( -csc( a_alpha(E) )^2 + cot( a_alpha(E) )/a_alpha(E) + 1/a_beta(E) );

x0 = input("Enter an initial guess: ");

% Run the algorithm
root = NRMethod(f_even, df_even, x0, 0.001, 0.001, bounds);

% Print the roots and plots
fprintf("The root: %f\n", root);

% Newton-Raphson Method Algorithm
function root = NRMethod(f, df, x0, tol_x, tol_y, bounds)
    
    x = x0;
    x_old = x - 2*tol_x;
    
    while abs(f(x)) > tol_y && abs(x - x_old) > tol_x
        
        x_old = x;
        
        if df(x) == 0.0
            x = x - tol_x;
        else
            x = x - f(x)/df(x);
        
            if ~bounds(x)
                x = NaN;
                break;
            end
        end
    
    end
    
    root = x;
end