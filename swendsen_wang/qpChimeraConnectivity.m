function [vertices, edges, A] = qpChimeraConnectivity(N,M,L)

Im = eye(M);
In = eye(N);
Il = eye(L);

e0110 = [ 0 1 ; 1 0 ];
e1000 = [ 1 0 ; 0 0 ];
e0001 = [ 0 0 ; 0 1 ];

Lm = eye(M) + diag(ones(1,M-1), 1) + diag(ones(1,M-1), -1);
Ln = eye(N) + diag(ones(1,N-1), 1) + diag(ones(1,N-1), -1);

A = kron(In, kron(Im, kron(e0110, ones(L)))) ...
  + kron(Lm, kron(In, kron(e1000, Il))) ...
  + kron(In, kron(Ln, kron(e0001, Il)));

vertices = (1:size(A, 1))';
[i, j] = ind2sub(size(A), find(triu(A) - eye(size(A))));
edges = cat(2, i, j);