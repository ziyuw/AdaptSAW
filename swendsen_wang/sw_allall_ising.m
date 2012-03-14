function [state, cprobs, energy_vec] = sw_allall_ising(state, ising_J, ising_h, iters, indices, edges, beta)
%SW_ALLALL_ISING Swendsen-Wang algorithm for fully-connect Ising model
%
%     [state, cprobs] = sw_allall_ising(state, ising_J, ising_h, iters)
%
% Inputs:
%        state DxN N D-dimensional {+1,-1} vectors
%      ising_J DxD This matrix should be symmetric, contributing energy: 
%                      - 0.5 \sum_i \sum_j J_{ij} s_i s_j
%                  That is both J_{ij} and J_{ji} contain the coupling between
%                  s_i and s_j. The 0.5 takes care of double-counting.
%      ising_h Dx1 or DxN
%                      field terms \sum_i h_i s_i
%                  There can be one common h, or one for each vector.
%        iters 1x1 Number of Markov chain steps to take.
%
% Outputs:
%        state DxN Markov chain state after an update that leaves the
%                  Ising model distribution(s) invariant.
%       cprobs DxN Marginal probabilities that each site would have been set to
%                  the "+1" state conditioned on the bond variables in the last update.

% Iain Murray, August 2007, November 2009

original = 0;
first_mod = 0;

potts_J = 2*ising_J;
if isvector(ising_h)
    ising_h = repmat(ising_h, 1, size(state, 2));
end

bond_probs = 1 - exp(-abs(potts_J));
flips = zeros(size(state));
flipprobs = zeros(size(state));

bonds = zeros(size(bond_probs));
if original || first_mod
    bonds_allowed = zeros(size(bond_probs));
end

energy_vec = zeros(iters,1);

for it = 1:iters
    energy_vec(it) = -state'*ising_J*state/(beta*2)-state'*ising_h/beta;

    if it/100 == ceil(it/100)
	fprintf('iteration: %3d in progress.\n', it);
    end

    for nn = 1:size(state,2)
        % Find bond variables
	if original
	    bonds_allowed = ((state(:,nn)*state(:,nn)').*potts_J > 0);
	    bonds = (rand(size(bond_probs)) < bond_probs) .* bonds_allowed;
	elseif first_mod
	    bonds_allowed(indices) = (state(edges(:,1),nn).*state(edges(:,2),nn).*potts_J(indices) > 0);
	    bonds(indices) = (rand(size(indices)) < bond_probs(indices));
	    bonds = bonds.* bonds_allowed;
	else
	    valid_entries = (state(edges(:,1),nn).*state(edges(:,2),nn).*potts_J(indices) > 0);
	    bonds(indices) = (rand(size(indices)) < bond_probs(indices));
	    bonds(indices) = bonds(indices).* valid_entries;
	end

        % Work out which spins to flip (and with what probability):
        flip_bias = 2*ising_h(:,nn).*state(:,nn); % bonus associated with sticking rather than flipping
        [flips(:,nn), flipprobs(:,nn)] = sw_allall_flip_conditionals(bonds, flip_bias);
    end
    flips = 2*flips - 1;

    if it == iters
        % This book-keeping has to be done before updating the state
        cprobs = flipprobs;
        downs = (state < 0);
        cprobs(downs) = 1 - cprobs(downs);
    end

    % New state
    state = state.*flips;
end
