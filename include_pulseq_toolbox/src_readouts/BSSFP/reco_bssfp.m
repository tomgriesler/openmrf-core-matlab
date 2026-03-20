%%
clear


cmaps = [];
[Images, PULSEQ, cmaps] = BSSFP_reco(['/home/ayde/University of Michigan Dropbox/Reina Ayde/data/pulseq_data/20260223/meas_MID00109_FID10093_260223_1311_reina_bssfp.dat'], cmaps);

xtv(Images);
