function [h, J, edges, indices, state] = param_setup(num_nodes, data, beta)

J = zeros(num_nodes);
for i =1:num_nodes
    edge_list = data.incident_list{i} + 1;
    nl = data.isg_edges(data.incident_list{i} + 1,:);
    for j = 1:size(nl, 1)
	J(nl(j, 1) + 1, nl(j, 2) + 1) = data.isg_J( edge_list(j) )*beta;
	J(nl(j, 2) + 1, nl(j, 1) + 1) = data.isg_J( edge_list(j) )*beta;
    end
end

edges = vertcat(data.isg_edges+1, data.isg_edges(:, [2,1])+1);
indices = (edges(:,1)-1)*num_nodes + edges(:,2);

h = data.isg_h*beta;
state = double(data.S_in);