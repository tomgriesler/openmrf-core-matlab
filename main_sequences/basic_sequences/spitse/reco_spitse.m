%%
clear
cmaps = [];
[Images, PULSEQ, study_info, cmaps, dcf2D, f0] = SPITSE_reco([], cmaps, [], []);
xtv(Images)
