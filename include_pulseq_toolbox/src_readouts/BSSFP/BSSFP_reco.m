function [Images, PULSEQ, cmaps] = BSSFP_reco(study, cmaps)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% ----- Input -----
% study:       enter path of study '.dat' or [] for dropdown select
% cmaps:       enter [] for calculating cmaps via openadapt or espirit
% zero_params: parameters for zero interpolation filling
% mod_reco:    0 -> sqrt(sum(abs(coils).^2,1)),
%              1 -> openadapt()
%              2 -> espirit()

% case: no input arguments
if nargin==0
    study       = [];
    cmaps       = [];
end



% read study info and pulseq workspace
[twix_obj, ~, PULSEQ] = pulseq_read_meas_siemens(study);
load('backup_260225_1608_workspace.mat');
%% import rawdata
% if PULSEQ.BSSFP.os_mode == 1
%     twix_obj.image.flagRemoveOS = 1;
% end
rawdata = squeeze(twix_obj.image());
[Nread, ~, ~] = size(rawdata);

%% reshape rawdata
N_iNAV = PULSEQ.BSSFP.N_iNAV;
N_rr   = PULSEQ.BSSFP.N_rr;
N_segment = PULSEQ.BSSFP.N_segment;
N_tot = N_iNAV + N_segment;
NCoils = size(rawdata,2);    
ksp_temp = zeros(Nread, NCoils, N_rr, N_segment);
inav_reorderd = zeros(Nread, NCoils, N_rr, N_iNAV);
for i = 1:N_rr
    for j = 1:N_tot
        if j > N_iNAV
            ksp_temp(:, :, i, j-N_iNAV) = rawdata(:, :, (N_tot*(i-1)) + j);
        elseif j <= N_iNAV
            inav_reorderd(:, :, i, j) = rawdata(:, :, (N_tot*(i-1)) + j);
        end
    end
end

%% reorder kspace data
traj = PULSEQ.BSSFP.caspr;
traj_names = fieldnames(traj);

ksp_reordered = zeros(Nread, NCoils, PULSEQ.FOV.Ny, PULSEQ.FOV.Nz);

for i = 1:length(traj_names)
    for j = 1:length(traj.(traj_names{i}))
        seg_coord = traj.(traj_names{i});
        ksp_reordered(:, :, seg_coord(j,1), seg_coord(j,2)) = ksp_temp(:, :, i, j);
    end
end



%% image reconstruction

% ifft
rawdata  = permute(ksp_reordered, [2, 1, 3, 4]);
images = fftshift(fft(fft(fft(rawdata, [], 4), [], 3), [], 2));
% calculate coil combined images
% cmaps = ones(size(images));

Images = squeeze(sqrt(sum(abs(images).^2,1)));

end

