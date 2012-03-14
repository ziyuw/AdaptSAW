function acf = acf_fft( x )

% Only the first few lags are "reliable!"
% This is a "true" autocorrelation, NOT a periodic/circular acf!
 
N = length( x );

% Make x into a column vector if it's not
x =   reshape( x, N, 1 );

mu = mean( x );

% Extend to pad with zeros
x_hat = [x-mu; zeros( N, 1) ];


X = fft( x_hat );
acf = ifft( X.*conj( X ) );

% Normalize by number of terms in the convolution:
acf = acf( 1: N );
acf = acf ./ ( N:-1:1 )';

% Scale by variance:
acf = acf  / acf(1);