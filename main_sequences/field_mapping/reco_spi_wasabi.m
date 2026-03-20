%% spiral reco
clear

zero_params.onoff  = 1;
zero_params.radius = 6.0;
zero_params.factor = 2.0;

cmaps = [];
[Images, PULSEQ, study_info, cmaps, dcf2D, f0] = SPI_reco([], cmaps, zero_params, [], []);
NImages = size(Images, 1);

%% mask and normalize data to S0
[mask_fit, mask_fit3D] = mg_get_mask_fit(squeeze(mean(abs(Images))), 'holes', NImages);
Images = Images .* mask_fit3D;

ImagesNorm = Images(1:end-1,:,:) * 0;
for j=1:PULSEQ.WASABI.n_ppm
    temp = squeeze(abs(Images(j,:,:)))./squeeze(abs(Images(end,:,:))); 
    ImagesNorm(j,:,:) = temp(:,:);
end 
clear temp;
ImagesNorm(isnan(ImagesNorm)) = 0;

xtv(ImagesNorm)

%% fit data:
[db1_Map, df0_Map, R2_Map, C_Map, D_Map] = WASABI_fit(ImagesNorm, mask_fit, PULSEQ.WASABI.f_off, PULSEQ.WASABI.tau, PULSEQ.WASABI.B1);

%% vis
WASABI_vis([-50 50], [0.75 1.25], 1, ImagesNorm, df0_Map, db1_Map, C_Map, D_Map, R2_Map, mask_fit, PULSEQ.WASABI.B1, PULSEQ.WASABI.f_off, PULSEQ.WASABI.tau);

%% vis interactive
WASABI_vis([-50 50], [0.75 1.25], 2, ImagesNorm, df0_Map, db1_Map, C_Map, D_Map, R2_Map, mask_fit, PULSEQ.WASABI.B1, PULSEQ.WASABI.f_off, PULSEQ.WASABI.tau);

