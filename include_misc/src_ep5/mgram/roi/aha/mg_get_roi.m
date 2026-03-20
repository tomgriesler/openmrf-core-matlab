function [roi] = mg_get_roi(I, settings, roi_old)

% Version: Maximilian Gram, 21.03.2024

repeat=1;

%% try old ROI
if ~isempty(roi_old)

    %% show map or reference image
    fig_no = 3456;
    figure(fig_no);
    set(gcf, 'Position', get(0, 'Screensize'));
    hold on
    imagesc(I, [settings.lb, settings.ub])
    colormap(settings.colormap);
    axis image;
    axis off;
    xlim(settings.xlims)    
    ylim(settings.ylims)    
    if isfield(settings, 'title')
        title(settings.title)
    end
    if isfield(settings, 'x_plot') && isfield(settings, 'y_plot')
        plot(settings.x_plot, settings.y_plot, 'r-', 'LineWidth', 1)
    end    
    
    %% repeat?
    plot(roi_old.xpos_smooth, roi_old.ypos_smooth, 'r-', 'LineWidth', 2)
    hold off
    answer = questdlg('repeat?', ...
        'repeat?', ...
        'repeat','next','next');
    switch answer
        case 'repeat'
            repeat = 1;
        case 'next'
            repeat = 0;
    end

    %% output
    if repeat==0
        roi.mask_smooth   = roi_old.mask_smooth;
        roi.mask_weighted = roi_old.mask_weighted;
        roi.xpos_raw      = roi_old.xpos_raw;
        roi.ypos_raw      = roi_old.ypos_raw;
        roi.xpos_smooth   = roi_old.xpos_smooth;
        roi.ypos_smooth   = roi_old.ypos_smooth;
    end

end

%% creat enew ROI
while repeat==1
clearvars -except I settings;

%% defaults
if ~isfield(settings, 'weight')
    settings.weight = 1;
end
if ~isfield(settings, 'lb')
   settings.lb = median(I,'all') - 3*std(std(I));
end
if ~isfield(settings, 'ub')
   settings.ub = median(I,'all') + 3*std(std(I));
end
if ~isfield(settings, 'xlims')
   settings.xlims = [1 size(I,2)];
end
if ~isfield(settings, 'ylims')
   settings.ylims = [1 size(I,1)];
end
if ~isfield(settings, 'colormap')
   settings.colormap = gray(1000);
end

%% show map or reference image
fig_no = 3456;
figure(fig_no);
set(gcf, 'Position', get(0, 'Screensize'));
hold on
imagesc(I, [settings.lb, settings.ub])
colormap(settings.colormap);
axis image;
axis off;
xlim(settings.xlims)    
ylim(settings.ylims)    
if isfield(settings, 'title')
    title(settings.title)
end
if isfield(settings, 'x_plot') && isfield(settings, 'y_plot')
    plot(settings.x_plot, settings.y_plot, 'r-', 'LineWidth', 1)
end

%% create roi using: imfreehand() or drawfreehand(), interp1(), smoothdata()

% get coordinates using drawfreehand()
myroi             = drawfreehand('FaceAlpha', 0, 'Multiclick', true);
mydlg             = warndlg('ROI completed?');
waitfor(mydlg);
xy_raw            = myroi.Position;

% periodic continuation
nxy                        = size(xy_raw,1);
xy_periodic                = zeros(2*nxy,2);
xy_periodic(1:nxy,1)       = xy_raw(:,1);
xy_periodic(1:nxy,2)       = xy_raw(:,2);
xy_periodic(nxy+1:2*nxy,1) = xy_raw(:,1);
xy_periodic(nxy+1:2*nxy,2) = xy_raw(:,2);

% smooth trajectory
filtersize     = round(nxy*0.20);
xy_smooth      = zeros(2*nxy,2);
xy_smooth(:,1) = smoothdata(xy_periodic(:,1),'sgolay',filtersize);
xy_smooth(:,2) = smoothdata(xy_periodic(:,2),'sgolay',filtersize);

% cut final trajectory
xy_final = xy_smooth(round(nxy/4):round(nxy/4)+nxy,:);

% output
xpos_raw    = xy_raw(:,1);
ypos_raw    = xy_raw(:,2);
xpos_smooth = xy_final(:,1);
ypos_smooth = xy_final(:,2);
mask_smooth = roipoly(I,xpos_smooth, ypos_smooth);

%% create weighted mask
if settings.weight == 1
    query_size       = 0.1;
    query_window_x   = floor(min(xpos_smooth))-0.5+query_size : query_size : ceil(max(xpos_smooth))+0.5;
    query_window_y   = floor(min(ypos_smooth))-0.5+query_size : query_size : ceil(max(ypos_smooth))+0.5;
    [meshx, meshy]   = meshgrid(query_window_x, query_window_y);
    query_x          = reshape(meshx, numel(meshx), 1);
    query_y          = reshape(meshy, numel(meshx), 1);
    [inside]         = inpolygon(query_x, query_y, xpos_smooth, ypos_smooth);
    query_map        = reshape(inside, size(meshx,1), size(meshx,2));
    query_map_mean   = zeros(size(query_map,1)/10,size(query_map,2)/10);
    for y = 1:size(query_map,1)/10
    for x = 1:size(query_map,2)/10
        x_start             = (x-1)*10+1;
        y_start             = (y-1)*10+1;
        temp(:,:)           = query_map( y_start:y_start+9 , x_start:x_start+9 );
        query_map_mean(y,x) = mean(mean(temp));
    end
    end
    mask_weighted = mask_smooth.*0;
    mask_weighted( floor(min(ypos_smooth)):ceil(max(ypos_smooth)), floor(min(xpos_smooth)):ceil(max(xpos_smooth)) ) = query_map_mean(:,:);
else
        mask_weighted = mask_smooth;
end
    
%% output
roi.mask_smooth   = mask_smooth;
roi.mask_weighted = mask_weighted;
roi.xpos_raw      = xpos_raw;
roi.ypos_raw      = ypos_raw;
roi.xpos_smooth   = xpos_smooth;
roi.ypos_smooth   = ypos_smooth;

%% repeat?
plot(xpos_smooth, ypos_smooth, 'r--', 'LineWidth', 2)
hold off
answer = questdlg('repeat?', ...
    'repeat?', ...
    'repeat','next','next');
switch answer
    case 'repeat'
        repeat = 1;
    case 'next'
        repeat = 0;
end
end

close(fig_no)

end

