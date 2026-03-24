function h = safe_plot(pns, gwr, dt, slice_orientation)
%SAFE_PLOT
%   Subplot 1–3: axis stimulation (X/Y/Z) + axis max + global ||p||
%   Subplot 4: Gx/Gy/Gz gradients
%
% Inputs:
%   pns, gwr : [N x 3]
%   dt       : [s]
%
% author: Maximilian Gram; University of Wuerzburg; 20.02.2026

if strcmp(slice_orientation, 'all')
    slice_orientation = 'axial  |  coronal  | sagittal';
end

N = size(pns,1);
t = (0:N-1) * dt * 1000; % ms

pnsnorm = sqrt(sum(pns.^2,2));
pGlob   = max(pnsnorm);
pAxis   = max(pns,[],1);

cXYZ = [0.0000 0.4470 0.7410;   % X
        0.8500 0.3250 0.0980;   % Y
        0.9290 0.6940 0.1250];  % Z

cTot  = [0 0 0];
lwSig = 2;
lwRef = 2;
lwTot = 1.5;

h = struct();
h.fig = figure('Color','w');
tl = tiledlayout(h.fig, 4, 1, 'TileSpacing','compact', 'Padding','compact');

dim = {'X','Y','Z'};

for k = 1:3
    ax = nexttile(tl,k); h.ax(k) = ax; %#ok<AGROW>
    hold(ax,'on'); grid(ax,'on');

    h.pns(k) = plot(ax, t, pns(:,k), '-', ...
        'Color', cXYZ(k,:), 'LineWidth', lwSig);

    h.maxAxis(k) = plot(ax, [t(1) t(end)], [pAxis(k) pAxis(k)], '--', ...
        'Color', cXYZ(k,:), 'LineWidth', lwRef);

    h.nrm(k) = plot(ax, t, pnsnorm, '--', ...
        'Color', cTot, 'LineWidth', lwTot);

    ylim(ax, [0 pGlob]);
    ylabel(ax, sprintf('%s stim [%%]', dim{k}));

    % remove x-axis labels & ticks to save space
    ax.XTickLabel = [];

    if k == 1        
        title(ax, ['Predicted PNS (global max ||p|| = ' num2str(round(pGlob,1)) '%  ->  ' slice_orientation]);
    end

    h.leg(k) = legend(ax, ...
        [h.pns(k) h.maxAxis(k) h.nrm(k)], ...
        {sprintf('%s stim', dim{k}), ...
         sprintf('%s max = %.0f%%', dim{k}, pAxis(k)), ...
         sprintf('||p|| (max %.0f%%)', pGlob)}, ...
        'Location','best');
end

ax4 = nexttile(tl,4); h.ax(4) = ax4;
hold(ax4,'on'); grid(ax4,'on');

h.gx = plot(ax4, t, gwr(:,1), '-', 'Color', cXYZ(1,:), 'LineWidth', lwSig);
h.gy = plot(ax4, t, gwr(:,2), '-', 'Color', cXYZ(2,:), 'LineWidth', lwSig);
h.gz = plot(ax4, t, gwr(:,3), '-', 'Color', cXYZ(3,:), 'LineWidth', lwSig);

ylabel(ax4, 'Gradients');
xlabel(ax4, 'Time [ms]');
h.leg(4) = legend(ax4, [h.gx h.gy h.gz], {'Gx','Gy','Gz'}, 'Location','best');

linkaxes(h.ax,'x');
xlim(h.ax, [t(1) t(end)]);
hZ = zoom(h.fig);
setAxesZoomMotion(hZ, h.ax(1), 'horizontal');

end