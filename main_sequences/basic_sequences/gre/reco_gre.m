%%
clear

zero_params.onoff  = 1;
zero_params.radius = 6.0;
zero_params.factor = 2.0;

cmaps = [];
[Images, PULSEQ, study_info, cmaps, f0] = GRE_reco([], cmaps, zero_params, []);

xtv(Images);
