function mu = roi_mean(map, xc, yc, r, r2)
 
    % define ROI with partial volume weighting
    N        = size(map,1);
    [xg, yg] = meshgrid(1:N, 1:N);    
    s        = 20;
    offs     = ((0:s-1)+0.5)/s - 0.5;
    [ox, oy] = meshgrid(offs, offs);
    ox       = ox(:)';
    oy       = oy(:)';    
    dx       = xg(:) - xc;
    dy       = yg(:) - yc;    
    in       = ((dx + ox).^2 + (dy + oy).^2) <= r^2;
    ROI      = reshape(mean(in, 2), N, N);
    ROI      = ROI / max(ROI(:));

    % filter low r2 values
    if nargin==5
        ROI(r2<0.95) = 0;
    end

    % calculate weighted mean
    vals    = map(ROI>0.05);
    weights = ROI(ROI>0.05);
    mu      = sum(vals .* weights) / sum(weights);

end