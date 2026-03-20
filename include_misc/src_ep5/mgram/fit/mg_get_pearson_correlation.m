function [correlation] = mg_get_pearson_correlation(x,y)

dx          = x - mean(x);
dy          = y - mean(y);
correlation = sum( dx .* dy ) / sqrt( sum( dx.^2 ) * sum( dy.^2 ) );

end

