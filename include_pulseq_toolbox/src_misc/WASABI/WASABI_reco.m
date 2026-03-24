function [db1_Map, df0_Map, R2_Map, C_Map, D_Map, ImagesNorm, mask_fit, PULSEQ, study_info] = WASABI_reco(path_raw, path_backup, vendor, cmaps, zero_params, mod_reco, ktraj_meas, mask_fit, flag_plot)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 23.03.2026

% ----- Input -----
% path_raw:    path of meas data or [] for select via uigetfile
% path_backup: path of pulseq workspace backup; not necessary for Siemens
% vendor:      vendor name
% cmaps:       enter [] for calculating cmaps via openadapt or espirit
% zero_params: parameters for zero interpolation filling
% mod_reco:    1 -> openadapt()
%              2 -> espirit()
% ktraj_meas:  [2 x Nid x Nadc] measured k-space trajectories
% mask_fit:    fit mask
% flag_plot:   0: no plot, 1: plot df0 and db1, 2: interactive fit results

% ----- Output -----
% db1_Map:    [ ]  fir result for B1+ deviation, 1 -> 100%
% df0_Map:    [Hz] fit result for off-resonance
% R2_Map:     R2 of fit
% C_map:      fit result C, see doi: 10.1002/mrm.26133
% D_Map:      fit result D, see doi: 10.1002/mrm.26133
% ImagesNorm: normalized WASABI images
% mask_fit:   fit mask
% PULSEQ:     struct with pulseq backup
% study_info: struct with header information

%% run spiral image reconstruction
[Images, PULSEQ, study_info] = SPI_reco(path_raw, path_backup, vendor, cmaps, zero_params, mod_reco, ktraj_meas);

%% mask and normalize data to S0
if isempty(mask_fit)
    mask_fit = mg_get_mask_fit(squeeze(mean(abs(Images))), 'holes');
end
ImagesNorm = Images(1:end-1,:,:) * 0;
for j=1:PULSEQ.WASABI.n_ppm
    temp = squeeze(abs(Images(j,:,:)))./squeeze(abs(Images(end,:,:))); 
    ImagesNorm(j,:,:) = temp .* mask_fit;
end 
clear temp;
ImagesNorm(isnan(ImagesNorm)) = 0;

%% fit data:
[db1_Map, df0_Map, R2_Map, C_Map, D_Map] = WASABI_fit(ImagesNorm, mask_fit, PULSEQ.WASABI.f_off, PULSEQ.WASABI.tau, PULSEQ.WASABI.B1);

%% vis
WASABI_vis([-50 50], [0.75 1.25], flag_plot, ImagesNorm, df0_Map, db1_Map, C_Map, D_Map, R2_Map, mask_fit, PULSEQ.WASABI.B1, PULSEQ.WASABI.f_off, PULSEQ.WASABI.tau);

end

