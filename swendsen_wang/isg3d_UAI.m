m = 9;
iters = 100000;
max_lag = 500;
load isg3d.mat;
if ~exist('file_loaded')
    load isg3d.mat;
    file_loaded = 1;
end

true_beta = beta_true;
true_beta = 1.0;
J = zeros(m^3);
for i =1:m^3
    edge_list = incident_list(i,:) + 1;
    nl = isg_edges(edge_list,:);
    for j = 1:size(nl, 1)
	J(nl(j, 1) + 1, nl(j, 2) + 1) = isg_J(edge_list(j))*true_beta;
	J(nl(j, 2) + 1, nl(j, 1) + 1) = isg_J(edge_list(j))*true_beta;
    end
end

[edges, indices, M] = isg3d_pattern(m);

h = isg_h;
state = double(S_in);

tic
[state, cprobs, energy_vec] = sw_allall_ising(state, J, h, iters, indices, edges, true_beta);
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
