function [Images, PULSEQ, study_info, cmaps, f0] = GRE_reco(study, cmaps, zero_params, mod_reco)

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
NImages = PULSEQ.GRE.n_TEs;
NRep    = PULSEQ.GRE.n_rep;
os_mode = PULSEQ.GRE.os_mode;

%% import rawdata
if os_mode == 1
    twix_obj.image.flagRemoveOS = 1;
end
rawdata = squeeze(twix_obj.image());
f0      = study_info.f0;

%% reshape rawdata
if NImages==1 && NRep==1
    Nx     = size(rawdata,1);
    NCoils = size(rawdata,2);    
    Ny     = size(rawdata,3);
    temp   = double(permute(rawdata, [2,3,1]));
    kSpace_raw(1,:,:,:) = temp(:,:,:);
    clear temp;
end

if NImages>1 || NRep>1
    Nx         = size(rawdata,1);
    NCoils     = size(rawdata,2);
    Ny         = size(rawdata,3)/NImages/NRep;
    kSpace_raw = zeros(NImages*NRep,NCoils,Ny,Nx);
    for nc = 1:NCoils
    for x  = 1:Nx
        temp = double(squeeze(rawdata(x,nc,:)));
        temp = reshape(temp,[Ny,NImages*NRep]);
        for ni = 1:NImages*NRep
            kSpace_raw(ni,nc,:,x) = squeeze(temp(:,ni));
        end
    end
    end
    clear temp nc x ni;
end

%% mean(image repetitions)
if NRep>1
    temp      = kSpace_raw;
    kSpace_raw = zeros(NImages,NCoils,Ny,Nx);
    for ni=1:NImages
        temp2(:,:,:,:) = temp( (ni-1)*NRep+1 : ni*NRep  ,:,:,:);
        temp2          = mean(temp2);
        kSpace_raw(ni,:,:,:) = temp2(1,:,:,:);
        clear temp2;
    end
    clear temp;
end

%% image reconstruction

% ifft
Images_coils = kspace2image(kSpace_raw);

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

