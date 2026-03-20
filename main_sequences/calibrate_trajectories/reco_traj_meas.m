%% get measured trajectory
clear

study_path = 'Q:/data/Pulseq/Rawdata/mgram/Skyra/250716_strukturphantom_git_test/';
study_name = 'meas_MID00736_FID228644_pulseq_traj_cmrf.dat';

[ktraj_meas, ktraj_hash] = TRAJ_reco([study_path study_name], [8:8:48]);
