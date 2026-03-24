%%
clear

zero_params.onoff  = 1;
zero_params.filter = 'kaiser';
zero_params.radius = 6.0;
zero_params.factor = 2.0;

[db1_Map, df0_Map, R2_Map, C_Map, D_Map, ImagesNorm, mask_fit, PULSEQ, study_info] = WASABI_reco([], [], 'Siemens', [], zero_params, 2, [], [], 1);

xtv(ImagesNorm)

%% vis
WASABI_vis([-50 50], [0.75 1.25], 1, ImagesNorm, df0_Map, db1_Map, C_Map, D_Map, R2_Map, mask_fit, PULSEQ.WASABI.B1, PULSEQ.WASABI.f_off, PULSEQ.WASABI.tau);

%% vis interactive
WASABI_vis([-50 50], [0.75 1.25], 2, ImagesNorm, df0_Map, db1_Map, C_Map, D_Map, R2_Map, mask_fit, PULSEQ.WASABI.B1, PULSEQ.WASABI.f_off, PULSEQ.WASABI.tau);

