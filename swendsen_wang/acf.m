function acf = compute_acf(seq, k)

n = size(seq,1);
mu = mean(seq);
va = var(seq);

acf = sum( (seq(1:(n-k)) - mu) .* (seq(k+1:n) - mu)) / ((n-k) * va);
