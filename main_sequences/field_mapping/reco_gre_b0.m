%%
clear

zero_params.onoff  = 1;
zero_params.radius = 6.0;
zero_params.factor = 2.0;

[Images, PULSEQ, study_info, cmaps, f0] = GRE_reco([], [], zero_params, []);
xtv(Images)

TEs = PULSEQ.GRE.TEs;
[df0_Map, dB0_Map, R2_Map, PhaseCorr] = mg_map_df0(Images, TEs, [], [-50 50], [1 0], 1);
