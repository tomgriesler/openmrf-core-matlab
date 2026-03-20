function [zq, xq, yq] = mg_get_heatmap(x_list, y_list, z_list, factor, show)

% Version: Maximilian Gram, 21.03.2024

if nargin<4
    factor = 1;
end

if nargin<5
    show = 0;
end

if size(x_list,1) < size(x_list,2)
    x_list = x_list';
end
if size(y_list,1) < size(y_list,2)
    y_list = y_list';
end
if size(z_list,1) < size(z_list,2)
    z_list = z_list';
end


xq_min = min(x_list);
xq_max = max(x_list);
yq_min = min(x_list);
yq_max = max(x_list);

Nq = floor(sqrt(size(x_list,1)) * factor);

[xq, yq] = meshgrid( linspace(xq_min,xq_max,Nq), linspace(yq_min,yq_max,Nq) );
zq       = griddata(x_list, y_list, z_list, xq, yq);

if show == 1
    figure(); imagesc(zq);
end

end