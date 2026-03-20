function [my_mean, my_std, my_median] = mg_get_weighted_mean(X, W)

% Version: Maximilian Gram, 21.03.2024

% create 1d array
x  = X(:);
w  = W(:);

% weighted mean
my_mean = sum( x.*w ) / sum(w);

% weighted standard deviation
my_std = sqrt( sum( w.*(x-my_mean).^2 ) / ( sum(w)-1 ) );

% median
w(w<1)    = 0;
x(w==0)   = [];
my_median = median(x);

end

