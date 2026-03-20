%% show measured kspace data
clear
[twix_obj, study_info, PULSEQ] = pulseq_read_meas_siemens();
rawdata = permute(twix_obj.image(), [2, 1, 3]);
xtv(rawdata)