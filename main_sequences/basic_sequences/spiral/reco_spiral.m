%% reco 2d
clear
cmaps = [];
[Images, PULSEQ, study_info, cmaps, dcf2D, f0, rawdata_noise] = SPI_reco([], cmaps, [], [], []);
xtv(Images)

%% reco 2d cg sense
cmaps = [];
[Images, PULSEQ, study_info, cmaps, dcf2D, f0] = SPI_reco_cg_sense([], cmaps, [], 25, [], []);
xtv(Images)

%% reco 3d stacked
clear
cmaps = [];
[Images, PULSEQ, study_info, cmaps, dcf2D, f0, rawdata_noise] = SPI_reco_3D_stacked([], cmaps, [], [], []);
xtv(Images)

%% reco 3d: only nearest neighbour approximation. feel free to improve :D
clear
cmaps = [];
[Images, PULSEQ] = SPI_reco_3D_NN();
xtv(Images)