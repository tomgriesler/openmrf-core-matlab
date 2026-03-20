function [xii, yii, ii] = get_poly_intersections(x1, y1, x2, y2)

n1 = numel(x1);
n2 = numel(x2);

ii = int16.empty(2,0)';
c  = 0;

for j=1:n1-1
for k=1:n2-1
    [xi, yi] = get_lin_intersection(x1(j), y1(j), x1(j+1), y1(j+1), x2(k), y2(k), x2(k+1), y2(k+1));  
    if ~isempty(xi)
        c      = c+1;
        xii(c) = xi;
        yii(c) = yi;
        ii_temp= [k j];
        ii = cat(1,ii, ii_temp); 
    end
end
end

if c==0
    xii = [];
    yii = [];
end

end

