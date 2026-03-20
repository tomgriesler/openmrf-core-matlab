%% show measured kspace data
clear
[twix_obj, study_info, PULSEQ] = pulseq_read_meas_siemens();
rawdata = permute(twix_obj.image(), [2, 1, 3]);
xtv(rawdata)

%% try stupid reco :D
images = kspace2image(mean(reshape(rawdata,size(rawdata,1),size(rawdata,2),PULSEQ.FOV.Ny,size(rawdata,3)/PULSEQ.FOV.Ny),4));
xtv(images)