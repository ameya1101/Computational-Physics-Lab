N_i = 10; %Initial dimension
N_f = 3000;
step = 5;

times = [];

dim = N_i;
while (dim <= N_f)
    
    % Generate the matrices A and B
    A = rand(dim, dim);
    B = rand(dim, 1);
    
    t_start = tic;
    X = A\B;
    t_stop = toc(t_start);
    
    times = [times, t_stop];
    
    dim = dim + step;
    
end

plot(N_i:step:N_f, times);
ylabel("Run times (s)", "interpreter", "latex");
xlabel("Dimension of the matrix $N$", "interpreter", "latex");