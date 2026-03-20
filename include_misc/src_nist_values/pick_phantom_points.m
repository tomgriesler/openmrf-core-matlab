function [xc, yc] = pick_phantom_points(N, ...
                                        map11, map12, map21, map22, ...
                                        clim11, clim12, clim21, clim22, ...
                                        cmap11, cmap12, cmap21, cmap22, ...
                                        xlims, ylims)
% pick_phantom_points
% Select N points with:
%   1) rough mouse click (ginput)
%   2) fine tuning via arrow keys (updates marker on ALL subplots)
%
% Display (2x2 grid):
%   (1,1) map11 with clim11 + cmap11
%   (1,2) map12 with clim12 + cmap12
%   (2,1) map21 with clim21 + cmap21
%   (2,2) map22 with clim22 + cmap22
%
% Workflow:
%   - Previously accepted points: blue M2
%   - Current point (fine-tuning): green M1 (moves on all subplots)
%
% Arrows        : move by 0.25
% Shift+Arrows  : move by 1
% Enter         : accept point
% Esc           : abort

xc = nan(N,1);
yc = nan(N,1);

M1 = 'rx';
M2 = 'gx';

for ii = 1:N

    % --- 1) plot maps (2x2) ---
    figure(333); clf;

    ax11 = subplot(2,2,1);
    imagesc(map11, clim11);
    axis image off;
    colormap(ax11, cmap11);
    xlim(xlims); ylim(ylims);
    title(sprintf('(1,1) Point %d/%d', ii, N), 'Interpreter','none');
    colorbar(ax11);

    ax12 = subplot(2,2,2);
    imagesc(map12, clim12);
    axis image off;
    colormap(ax12, cmap12);
    xlim(xlims); ylim(ylims);
    title(sprintf('(1,2) Point %d/%d', ii, N), 'Interpreter','none');
    colorbar(ax12);

    ax21 = subplot(2,2,3);
    imagesc(map21, clim21);
    axis image off;
    colormap(ax21, cmap21);
    xlim(xlims); ylim(ylims);
    title(sprintf('(2,1) Point %d/%d', ii, N), 'Interpreter','none');
    colorbar(ax21);

    ax22 = subplot(2,2,4);
    imagesc(map22, clim22);
    axis image off;
    colormap(ax22, cmap22);
    xlim(xlims); ylim(ylims);
    title(sprintf('(2,2) Point %d/%d', ii, N), 'Interpreter','none');
    colorbar(ax22);

    % --- 1b) plot previously accepted points (blue) on all subplots ---
    if ii > 1
        xp = xc(1:ii-1);
        yp = yc(1:ii-1);

        axes(ax11); hold on; plot(xp, yp, M2, 'MarkerSize', 12, 'LineWidth', 2);
        axes(ax12); hold on; plot(xp, yp, M2, 'MarkerSize', 12, 'LineWidth', 2);
        axes(ax21); hold on; plot(xp, yp, M2, 'MarkerSize', 12, 'LineWidth', 2);
        axes(ax22); hold on; plot(xp, yp, M2, 'MarkerSize', 12, 'LineWidth', 2);
    else
        % ensure hold state consistent for current marker
        axes(ax11); hold on;
        axes(ax12); hold on;
        axes(ax21); hold on;
        axes(ax22); hold on;
    end

    % --- 2) rough position (click in ax12 by default) ---
    axes(ax12);
    [x, y] = ginput(1);

    % clamp to limits immediately
    x = min(max(x, xlims(1)), xlims(2));
    y = min(max(y, ylims(1)), ylims(2));

    % --- 2b) current point marker (green) on all subplots ---
    axes(ax11); h11 = plot(x, y, M1, 'MarkerSize', 15, 'LineWidth', 2);
    axes(ax12); h12 = plot(x, y, M1, 'MarkerSize', 15, 'LineWidth', 2);
    axes(ax21); h21 = plot(x, y, M1, 'MarkerSize', 15, 'LineWidth', 2);
    axes(ax22); h22 = plot(x, y, M1, 'MarkerSize', 15, 'LineWidth', 2);

    title(ax12, ...
        'Arrows = 0.25 | Shift = 1 px | Enter = accept | Esc = abort', ...
        'Interpreter','none');
    drawnow;

    % --- 3) fine tuning loop ---
    while true
        waitforbuttonpress;
        key  = get(gcf, 'CurrentKey');
        mods = get(gcf, 'CurrentModifier');

        step = 0.25;
        if any(strcmp(mods, 'shift'))
            step = 1;
        end

        if strcmp(key, 'return')
            break;
        elseif strcmp(key, 'escape')
            xc = xc(1:ii-1);
            yc = yc(1:ii-1);
            close(333);
            warning('Aborted after %d points.', numel(xc));
            return;
        elseif strcmp(key, 'leftarrow')
            x = x - step;
        elseif strcmp(key, 'rightarrow')
            x = x + step;
        elseif strcmp(key, 'uparrow')
            y = y - step;   % imagesc coordinate system
        elseif strcmp(key, 'downarrow')
            y = y + step;
        else
            continue;
        end

        % clamp to limits
        x = min(max(x, xlims(1)), xlims(2));
        y = min(max(y, ylims(1)), ylims(2));

        % update current marker on ALL subplots
        set(h11, 'XData', x, 'YData', y);
        set(h12, 'XData', x, 'YData', y);
        set(h21, 'XData', x, 'YData', y);
        set(h22, 'XData', x, 'YData', y);
        drawnow limitrate;
    end

    % store accepted point
    xc(ii) = x;
    yc(ii) = y;
end

close(333);

% --- 4) final plot: show all 2x2 with accepted points (blue) ---
figure(334);

ax11 = subplot(2,2,1);
imagesc(map11, clim11); axis image off;
colormap(ax11, cmap11);
xlim(xlims); ylim(ylims);
title('(1,1) Final', 'Interpreter','none');
colorbar(ax11); hold on;
plot(xc, yc, M2, 'MarkerSize', 12, 'LineWidth', 2);

ax12 = subplot(2,2,2);
imagesc(map12, clim12); axis image off;
colormap(ax12, cmap12);
xlim(xlims); ylim(ylims);
title('(1,2) Final', 'Interpreter','none');
colorbar(ax12); hold on;
plot(xc, yc, M2, 'MarkerSize', 12, 'LineWidth', 2);

ax21 = subplot(2,2,3);
imagesc(map21, clim21); axis image off;
colormap(ax21, cmap21);
xlim(xlims); ylim(ylims);
title('(2,1) Final', 'Interpreter','none');
colorbar(ax21); hold on;
plot(xc, yc, M2, 'MarkerSize', 12, 'LineWidth', 2);

ax22 = subplot(2,2,4);
imagesc(map22, clim22); axis image off;
colormap(ax22, cmap22);
xlim(xlims); ylim(ylims);
title('(2,2) Final', 'Interpreter','none');
colorbar(ax22); hold on;
plot(xc, yc, M2, 'MarkerSize', 12, 'LineWidth', 2);

sgtitle('Final selected points', 'Interpreter','none');

end