function [edges, indeces, M] = isg3d_pattern(m)

num_edges = 6*m^3;
edges = zeros(num_edges, 2);

counter = 1;
for i = 1:m^3
    l = ceil(i/(m^2));
    f = ceil((i-(l-1)*m^2)/m);

    if (i-(l-1)*m^2)/m == f
	edges((i-1)*6+1, :) = [i, i-m+1];
    else
	edges((i-1)*6+1, :) = [i, i+1];
    end

    if i == (l-1)*m^2+(f-1)*m+1
	edges((i-1)*6+2, :) = [i, i+m-1];
    else
	edges((i-1)*6+2, :) = [i, i-1];
    end

    if f == 1
	edges((i-1)*6+3, :) = [i, i+(m-1)*m];
    else
	edges((i-1)*6+3, :) = [i, i-m];
    end

    if f == m
	edges((i-1)*6+4, :) = [i, i-(m-1)*m];
    else
	edges((i-1)*6+4, :) = [i, i+m];
    end

    if l == 1
	edges((i-1)*6+5, :) = [i, i+(m-1)*m^2];
    else
	edges((i-1)*6+5, :) = [i, i-m^2];
    end

    if l == m
	edges((i-1)*6+6, :) = [i, i-(m-1)*m^2];
    else
	edges((i-1)*6+6, :) = [i, i+m^2];
    end
end

indeces = (edges(:,1)-1)*m^3 + edges(:,2);

M = zeros(m^3);
M(indeces) = ones(num_edges, 1);
