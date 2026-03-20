function [roi] = mg_get_roi_lv(I, settings, roi_old)

% Version: Maximilian Gram, 21.03.2024

% ----- input: -----
% I:           relaxation time map or magnitude image
% settings:    structured object for plot settings
%  .lb:        lower boundary for colorbar
%  .ub:        upper boundary for colorbar
%  .colormap:  rgb colormap
%  .xlims:     limits for x axis
%  .ylims:     limits for y axis
% roi_old:     old rois for visualization
%  .endo:      endocardial roi
%  .epi:       epicardial  roi

%% load endocardial and epicardial rois
if nargin<3
    roi_old = [];
end
if isempty(roi_old)
    roi_old_endo = [];
    roi_old_epi  = [];
else
    roi_old_endo = roi_old.endo;
    roi_old_epi  = roi_old.epi;
end

%% create endo/epicardial rois and masks

% epicardial roi
settings.title = 'create ROI for epicardium (out)';
roi.epi        = mg_get_roi(I, settings, roi_old_epi);

% endocardial roi
settings.title  = 'create ROI for endocardium (blood pool)';
settings.x_plot = roi.epi.xpos_smooth;
settings.y_plot = roi.epi.ypos_smooth;
roi.endo        = mg_get_roi(I, settings, roi_old_endo);
    
% left ventricular roi
roi.mask_smooth   = roi.epi.mask_smooth   - roi.endo.mask_smooth;
roi.mask_weighted = roi.epi.mask_weighted - roi.endo.mask_weighted;

% calculate myocard center and angles of septal borders
[sy, sx]           = size(I);
[meshx,meshy]      = meshgrid(1:sx, 1:sy);
roi.yc             = sum(sum(meshy.*roi.mask_weighted)) / sum(sum(roi.mask_weighted));
roi.xc             = sum(sum(meshx.*roi.mask_weighted)) / sum(sum(roi.mask_weighted));
clear sx sy meshx meshy;

end



