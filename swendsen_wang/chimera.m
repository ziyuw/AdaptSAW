iters = 100000;
max_lag = min(500, iters-1);

data = load('Chimera.mat');

beta = 1.0;
num_nodes = 128;

[h, J, edges, indices, state] = param_setup(num_nodes, data, beta);


tic	
[state, cprobs, energy_vec] = sw_allall_ising(state, J, h, iters, indices, edges, beta);
toc

figure;
plot(energy_vec);
grid on;

figure;
as = zeros(max_lag, 1);
for i = 0:max_lag;
    as(i+1) = acf(energy_vec, i);
end
plot(as);
axis([0 max_lag 0.0 1.0])
grid on;
