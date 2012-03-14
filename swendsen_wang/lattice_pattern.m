function [edges, indeces, M] = lattice_pattern(m)

upper = m^2;
num_edges = ((m*(m-1))*2 + m*2)*2;
edges = zeros(num_edges, 2);

counter = 1;
for i = 1:m^2
    f = ceil(i/m);
    if i/m == f
	edges((i-1)*4+1, :) = [i, i-m+1];
    else
	edges((i-1)*4+1, :) = [i, i+1];
    end

    if i == (f-1)*m+1
	edges((i-1)*4+2, :) = [i, i+m-1];
    else
	edges((i-1)*4+2, :) = [i, i-1];
    end

    if f == 1
	edges((i-1)*4+3, :) = [i, i+(m-1)*m];
    else
	edges((i-1)*4+3, :) = [i, i-m];
    end

    if f == m
	edges((i-1)*4+4, :) = [i, i-(m-1)*m];
    else
	edges((i-1)*4+4, :) = [i, i+m];
    end
end

indeces = (edges(:,1)-1)*m^2 + edges(:,2);

M = zeros(m^2);
M(indeces) = ones(num_edges, 1);
