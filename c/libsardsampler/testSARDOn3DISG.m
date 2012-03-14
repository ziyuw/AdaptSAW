% NOTE RUNS SARDONICS2 IN V9 DIRECTORY

beep off;

% FOR DEBUG
%dbmex on;

% First, set all paths by running addPaths.m in enforcedTypicality
curr_dir = pwd;
cd ~/project/bpstuff/snapSampling/enforcedTypicality
addPaths;
addpath('~/project/bpstuff/snapSampling/SARDONICS2/ccode/V9');

cd (curr_dir);

L = 10;
M = 10;
N = 10;

disp('Making Graph...');
[edges_i, nbrs_list_i, incident_edge_list_i] = makeCubeGraphPeriodic_C(L,M,N);
disp('Done');

num_nodes = L*M*N;
num_edges = size( edges_i, 1);

load CRAP;
%load GAUSS_CRAP;
%J = makeRandomJEA( num_edges );
h = zeros( num_nodes,1);

% Freezing beta:
beta_crit = 0.9;

% beta_true = beta_crit;
beta_true = 1.2;

%gamma_lo = beta_true;
%gamma_hi = 0.89;

gamma_lo = beta_true;
gamma_hi = 0.89;

P_LH = 0.005;
P_HL = 0.04;
%P_LH = 0.0;
%P_HL = 0.0;
P_LL = 1-(P_LH+P_HL);

% The initial state is random:
S_in = ones( num_nodes, 1 );
S_in( find(rand(num_nodes,1) < 0.5) ) = -1;

% % 16.Dec.10: KEEP
% iters_per_move=70;
% SAWlength_min = 1;
% SAWlength_max = 1;

iters_per_move=70;
SAWlength_min = 1;
SAWlength_max = 1;

evaluateOverlap = 1;

%n_moves = 2e4;
n_moves = 1e3;

[h_effective_in] = getEffectiveFields_C(edges_i, incident_edge_list_i, J, h, S_in);

% Initial energy:
E_in = -double(S_in')*h - double(S_in')*(h_effective_in-h)/2; 

n_trials = 1;
burn_in = 1000;
acf_sard = zeros( n_moves-burn_in+1, 1 );

tic;
  
  for trial=1:n_trials
   
    disp(['On trial ' num2str(trial) ' of ' num2str( n_trials )]);
    
    % For SARDONICS:
    [E_samples_vectors, E_proposed_vectors, log_f_fwd_vectors, log_f_rev_vectors, MH_ratio_vectors, move_type_vectors, overlap_vector] = SARDRun( nbrs_list_i, incident_edge_list_i, edges_i, J, h, int32(S_in), h_effective_in, beta_true, gamma_lo, gamma_hi, int32(SAWlength_min), int32(SAWlength_max), int32(iters_per_move), int32(n_moves), P_LL, P_LH, P_HL, int32(evaluateOverlap) );
   
    if evaluateOverlap == 1 
      acf_sard = acf_sard + acf_fft( double(overlap_vector(burn_in:end)) );
    end
      
  end
  acf_sard = acf_sard / n_trials;
  
toc;

% Plot the first replica's delta E over move types
figure;
% LL-moves:
I = find( move_type_vectors(:,1) == 0 );
plot( I, E_proposed_vectors(I,1)-E_samples_vectors(I,1) );

figure;
% LH-moves:
I = find( move_type_vectors(:,1) == 1 );
plot( I, E_proposed_vectors(I,1)-E_samples_vectors(I,1) );

figure;
% HL-moves:
I = find( move_type_vectors(:,1) == 2 );
plot( I, E_proposed_vectors(I,1)-E_samples_vectors(I,1) );

% scatter( -beta_true*(E_proposed_vectors(I,1)-E_samples_vectors(I,1)), log_f_rev_vectors(I,1)-log_f_fwd_vectors(I,1) );

% II = find( E_proposed_vectors(I,1)-E_samples_vectors(I,1) <= 0 );
% plot( sort( MH_ratio_vectors(I(II),1)) );