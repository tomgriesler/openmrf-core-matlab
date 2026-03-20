function [xi, yi] = get_lin_intersection(x11, y11, x12, y12, x21, y21, x22, y22)

%% no common range
if ( max(x11,x12)<min(x21,x22) || max(x21,x22)<min(x11,x12) || max(y11,y12)<min(y21,y22) || max(y21,y22)<min(y11,y12)  )
    xi = [];
    yi = [];
%% search for intersection of two linear functions
else
% lin fun 1: (x11, y11) --- (x12, y12)
% lin fun 2: (x21, y21) --- (x22, y22)

m1 = (y12-y11) / (x12-x11);
m2 = (y22-y21) / (x22-x21);
t1 = y11 - m1 * x11;
t2 = y21 - m2 * x21;

if (m1~=m2) % not parallel
    xi = (t2-t1) / (m1-m2);
    yi = m1 * xi + t1;
else % parallel
    xi = [];
    yi = [];
end

% intersection outside line range
xrange = [x11 x12 x21 x22];
xrange = sort(xrange);
xrange = [ xrange(2) xrange(3) ];
yrange = [y11 y12 y21 y22];
yrange = sort(yrange);
yrange = [ yrange(2) yrange(3) ];

if ( xi<xrange(1) || xrange(2)<xi || yi<yrange(1) || yrange(2)<yi )
    xi = [];
    yi = [];
end    
    
end

end

