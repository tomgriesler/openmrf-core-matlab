%% show measured kspace data
clear
[twix_obj, study_info, PULSEQ] = pulseq_read_meas_siemens();
rawdata = permute(twix_obj.image(), [2, 3, 1]);
xtv(rawdata)

%% try reco
NCoils  = size(rawdata, 1);
ktraj   = PULSEQ.ktraj_reco;
ktraj   = ktraj(:,:);
rawdata = rawdata(:,:);
reco    = mg_mixed_mrf_reco(rawdata, ktraj, PULSEQ.FOV.Nxy, PULSEQ.FOV.fov_xy);