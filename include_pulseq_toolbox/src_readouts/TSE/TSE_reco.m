function [Images, PULSEQ, study_info, cmaps, f0] = TSE_reco(study, cmaps, zero_params, mod_reco)

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
    mod_reco = 2;
end

% read study info and pulseq workspace
[twix_obj, study_info, PULSEQ] = pulseq_read_meas_siemens(study);

%% read parameters for reconstruction
NImages = 1;
NTurbo  = PULSEQ.TSE.n_echo;
os_mode = PULSEQ.TSE.os_mode;

%% import/sort kSpace rawdata
if os_mode == 1
    twix_obj.image.flagRemoveOS = 1;
end
rawdata = squeeze(twix_obj.image());
f0      = study_info.f0;

%% reshape rawdata
rawdata    = double(permute(rawdata, [2,3,1]));
Nx         = size(rawdata,3);
Ny         = size(rawdata,2)/NImages;
NCoils     = size(rawdata,1);
kSpace_raw = zeros(NCoils, NImages, Ny, Nx);

for j=1:NCoils
    temp(:,:) = rawdata(j,:,:);
    k_start   = 1;
    k_end     = Ny;
    for k=1:NImages
        kSpace_raw(j,k,:,:) = rawdata(j,k_start:k_end,:);
        k_start = k_start + Ny;
        k_end   = k_end + Ny;        
    end
end
clear temp k_start k_end

PEorder  = PULSEQ.TSE.PEorder;
PEorder  = PEorder + Ny/2 + 1;
PEorder  = reshape(PEorder,[Ny,1]);

kSpace_sorted = zeros(NCoils,NImages,Ny,Nx);
for j=1:NCoils
for k=1:NImages
    temp(:,:) = kSpace_raw(j,k,:,:);
    for l=1:Ny
        kSpace_sorted(j,k,PEorder(l),:) = temp(l,:);        
    end    
end
end

%% image reconstruction

% ifft
Images_coils = kspace2image(permute(kSpace_sorted, [2,1,3,4]));

% calculate cmaps
if isempty(cmaps)
    if NImages==1
        temp = squeeze(Images_coils);
    else
        temp = squeeze(mean(Images_coils));
    end
    if mod_reco==1
        [~, cmaps] = openadapt(temp);    
    end
    if mod_reco==2
        cmaps = mg_espirit_cmaps(temp, 0.02, 0.95, 24, [6,6]);
    end
    clear temp;
end

% calculate coil combined images
Images = zeros(NImages, Ny, Nx);
for j=1:NImages
    if mod_reco == 0
        temp          = squeeze(Images_coils(j,:,:,:));
        temp          = squeeze(sqrt(sum(abs(temp).^2,1)));
        Images(j,:,:) = temp(:,:);
        clear temp;
    end    
    if mod_reco>1   
        Images(j,:,:) = squeeze(sum( squeeze(Images_coils(j,:,:,:)) .* conj(cmaps) ));
    end    
end

%% zero filling
if zero_params.onoff == 1
    Images = mg_zero_filling(Images, zero_params);
end

end

