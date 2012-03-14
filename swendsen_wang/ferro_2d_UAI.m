m = 60;
iters = 100000;
max_lag = 500;
true_beta = 1/2.27;
true_beta = 1;
[edges, indices, M] = lattice_pattern(m);

J = ((ones(m^2)-eye(m^2)).*M)*(true_beta);
h = zeros(m^2, 1);

state = ones(m^2,1);
state(1:(m^2)/2,1) = state(1:(m^2)/2,1)*(-1);
state = state(randperm(m^2));
figure(1)
tic
[state, cprobs, energy_vec] = sw_allall_ising(state, J, h, iters, indices, edges, true_beta);
toc
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

