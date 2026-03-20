function [mask] = mg_get_aha_segm_mask(I, yc, xc, ac, ad)

% Version: Maximilian Gram, 21.03.2024

% xc: x coord center
% yc: y coord center
% ac: angle center
% ad: angle difference -> segment size

[sy, sx] = size(I);
mask     = zeros(sy,sx);

for y=1:sy
for x=1:sx
    alpha = wrapTo2Pi(atan2(y-yc, x-xc));
    diff = min( mod(ac-alpha,2*pi), mod(alpha-ac,2*pi) );
    if diff <= ad
        mask(y,x) = 1;
    end
end
end


end

