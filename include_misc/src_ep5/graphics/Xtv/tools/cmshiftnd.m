function x = cmshiftnd( x, shifts)
%Function to circularly shift N-D arrays
%M.A. Griswold 9/11/98

% minor tweaks:
%   *   x = f(x) may be faster / more RAM-economic
%   *   return immediately if nargin < 2 or shifts is all zeros
% M. Völker 03.12.2012

    if nargin < 2 || all(shifts(:) == 0)
       return                       % no shift
    end

    sz      = size( x );
    numDims = ndims(x);             % number of dimensions
    idx = cell(1, numDims);         % creates cell array of empty matrices,
                                    % one cell for each dimension

    for k = 1:numDims               %#ok <--- no parfor-hint

        m = sz(k);
        p = ceil(shifts(k));

        if p < 0
            p = m + p;
        end

        idx{k} = [p+1:m  1:p];
    end

    % Use comma-separated list syntax for N-D indexing.
    x = x(idx{:});

end
