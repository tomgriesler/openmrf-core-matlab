function [Images, PULSEQ, cmaps] = BSSFP_reco(path_raw, path_backup, vendor, cmaps, zero_params, mod_reco)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V2, 22.03.2026; unify for different vendors

% ----- Input -----
% path_raw:    path of meas data or [] for select via uigetfile
% path_backup: path of pulseq workspace backup; not necessary for Siemens
% vendor:      vendor name
% cmaps:       enter [] for calculating cmaps via openadapt or espirit
% zero_params: parameters for zero interpolation filling
% mod_reco:    0 -> sqrt(sum(abs(coils).^2,1)),
%              1 -> openadapt()

% case: no input arguments
if nargin==0
    path_raw    = [];
    path_backup = [];
    vendor      = [];
    cmaps       = [];
    zero_params = [];
    mod_reco    = [];
end

% defaults
if isempty(zero_params)
    zero_params.onoff  = 1;
    zero_params.radius = 6.0;
    zero_params.factor = 2.0;
end
if isempty(mod_reco)
    mod_reco = 1;
end

% import rawdata, read study info and load pulseq workspace
[rawdata, ~, PULSEQ, study_info] = pulseq_read_meas(path_raw, path_backup, vendor);

%% read dimensions
[NCoils, ~, Nx] = size(rawdata);
Ny = PULSEQ.FOV.Ny;
Nz = PULSEQ.FOV.Nz;

%% sort k-space
kspace = zeros(Nz, NCoils, Ny, Nx);
for nz = 1:Nz
    temp = rawdata(:,(nz-1)*Ny+1:nz*Ny,:);
    kspace(nz,:,:,:) = temp;
end
clear temp;
kspace = permute(kspace, [2,1,3,4]);

%% ifft reconstruction and coil combination
images_coils = fftshift(fft(fft(fft(kspace, [], 4), [], 3), [], 2));

if mod_reco==0
    Images = squeeze(sqrt(sum(abs(images_coils).^2,1)));
    cmaps  = [];
else
    [Images, cmaps] = openadapt(images_coils);
end

%% zero filling
if zero_params.onoff == 1
    Images = mg_zero_filling(Images, zero_params);
end

end

