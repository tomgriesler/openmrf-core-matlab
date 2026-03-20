function [Images, PULSEQ, study_info, cmaps, dcf2D, f0, rawdata_noise] = SPI_reco_3D(study, cmaps, zero_params, mod_reco, ktraj_meas)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% ----- Input -----
% study:       enter path of study '.dat' or [] for dropdown select
% cmaps:       enter [] for calculating cmaps via openadapt or espirit
% zero_params: parameters for zero interpolation filling
% mod_reco:    1 -> openadapt()
%              2 -> espirit()
% ktraj_meas:  measured k-space trajectories

% case: no input arguments
if nargin==0
    study       = [];
    cmaps       = [];
    zero_params = [];
    mod_reco    = [];
    ktraj_meas  = [];
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

%% import spiral rawdata and noise pre-scans
if isfield(PULSEQ.SPI, 'Nnoise')
    Nnoise = PULSEQ.SPI.Nnoise;  
else
    Nnoise = 0;
end
[rawdata1, NImages, NRead, NCoils, rawdata_noise] = SPI_get_rawdata(twix_obj, Nnoise);
NR            = PULSEQ.SPI.NR;
Nz            = PULSEQ.FOV.Nz;
NImages       = NImages / NR / Nz;
f0            = study_info.f0;
rawdata1      = permute(rawdata1, [3, 1, 2]);
rawdata_noise = permute(rawdata_noise, [3, 1, 2]);

%% noise pre-whitening
if ~isempty(rawdata_noise)
    [rawdata1] = mg_noise_prewhitening(rawdata1, rawdata_noise, 'cholesky', 1);
end

%% import measured k-space trajectories
if isempty(ktraj_meas)
    ktraj_reco = SPI_load_ktraj(PULSEQ);
    if size(ktraj_reco,1) == 3
        ktraj_reco(3,:,:) = [];
    end
else
    ktraj_reco = ktraj_meas;
end
if isfield(PULSEQ.SPI, 'proj')
    ktraj_reco = ktraj_reco(:,PULSEQ.SPI.proj.id,:);    
elseif isfield(PULSEQ.SPI, 'phi_id')
    ktraj_reco = ktraj_reco(:,PULSEQ.SPI.phi_id,:);    
end

%% adc padding
if isfield(PULSEQ.SPI, 'adcNPad')
    if numel(PULSEQ.SPI.adcNPad) == 2
        rawdata1   = rawdata1(:,:, PULSEQ.SPI.adcNPad(1):PULSEQ.SPI.adcNPad(2));
        ktraj_reco = ktraj_reco(:,:, PULSEQ.SPI.adcNPad(1):PULSEQ.SPI.adcNPad(2));
    else
        rawdata1   = rawdata1(:,:,PULSEQ.SPI.adcNPad+1:end);
        ktraj_reco = ktraj_reco(:,:,PULSEQ.SPI.adcNPad+1:end);
    end
    NRead = size(rawdata1, 3);
end

%% separate NRs

if NImages==1
    rawdata2(:,1,:,:) = rawdata1(:,:,:);
else
    rawdata2 = zeros(NCoils, NImages, Nz*NR, NRead);
    for ni=1:NImages
        rawdata2(:,ni,:,:) = rawdata1(:,(ni-1)*Nz*NR+1 : ni*Nz*NR,:);
    end
end
clear rawdata1;

rawdata3 = zeros(NCoils, NImages, Nz, NR, NRead);
for nc = 1:NCoils
for ni = 1:NImages
for nz = 1:Nz
    temp = squeeze(rawdata2(nc, ni, (nz-1)*NR+1 : nz*NR, : ));
    rawdata3(nc, ni, nz, :, :) = temp(:,:);
end
end
end
clear nc ni nz temp rawdata2;

%% ifft for kz direction
rawdata4 = zeros(NCoils, NImages, Nz, NR, NRead);
for nc = 1:NCoils
for ni = 1:NImages
    temp = squeeze(rawdata3(nc,ni,:,:,:));
    temp = kspace2image(temp, [1, 0, 0]);
    rawdata4(nc,ni,:,:,:) = temp(:,:,:);
end
end
clear nc ni temp rawdata3;

%% join NRs
rawdata5 = zeros(Nz, NCoils, NImages, NRead*NR);
for j  = 1:NCoils
for k  = 1:NImages
for nz = 1:Nz    
    temp = squeeze(rawdata4(j,k,nz,:,:));
    temp = temp(:);
    rawdata5(nz, j, k, :) = temp(:);
end
end
end
clear j k nz temp rawdata4;

%% set nufft prams for reconstruction
Nxy = PULSEQ.FOV.Nxy;
N   = [Nxy, Nxy];         % [ ]   matrix size
fov = PULSEQ.FOV.fov_xy;  % [m]   fov size
kx  = ktraj_reco(1,:,:);  % [1/m] x-trajecotry
ky  = ktraj_reco(2,:,:);  % [1/m] y-trajecotry
kx  = kx(:);              % join NRs
ky  = ky(:);              % join NRs
kxy = [kx, ky];           % Nadc x 2

%% calculate density compensation function: DCF
nufft_args = {N, [6 6], 2*N, N/2, 'table', 2^12, 'minmax:kb'};
nufft_st   = nufft_init( kxy / max(abs(kxy(:)))*pi, N, [6 6], N*2, N/2, 'minmax:kb');
G          = Gmri(kxy, true(N), 'fov', fov, 'basis', {'rect'}, 'nufft', nufft_args);
dcf        = abs(mri_density_comp( kxy, 'pipe', 'G', G.arg.Gnufft)) *2^2 *Nxy^2;
dcf2D      = mg_get_heatmap( kx, ky, dcf, 2 );

%% calculate nonuniform fast fourier transform: NUFFT
images_coils = zeros(Nz, NImages, NCoils, Nxy, Nxy);
parfor nz = 1:Nz
   for ni = 1:NImages
    temp_k                    = squeeze(rawdata5(nz,:,ni,:));
    temp_reco                 = permute( nufft_adj( temp_k.' .* repmat(dcf, [1, NCoils]), nufft_st ), [3, 1, 2]);
    images_coils(nz,ni,:,:,:) = temp_reco(:,:,:);
end
end
clear rawdata5;

%% calculate coil sensitivity maps: cmaps
if isempty(cmaps)
    cmaps = zeros(Nz, NCoils, Nxy, Nxy);
    for nz=1:Nz    
        if mod_reco==1
            if NImages>1
                [~, temp_cmaps] = openadapt(squeeze(mean(images_coils(nz,:,:,:,:))));
            else
                [~, temp_cmaps] = openadapt(squeeze(images_coils(nz,:,:,:,:)));
            end
        end
        if mod_reco==2
            if NImages>1
                temp_cmaps = mg_espirit_cmaps(squeeze(mean(images_coils(nz,:,:,:,:))), 0.02, 0.95, 24, [6,6]);
            else
                temp_cmaps = mg_espirit_cmaps(squeeze(images_coils(nz,:,:,:,:)), 0.02, 0.95, 24, [6,6]);
            end
        end
        cmaps(nz,:,:,:) = temp_cmaps(:,:,:);
    end
end
clear temp_cmaps nz;

%% calculate coil combined images
Images_raw = zeros(NImages, Nz, Nxy, Nxy);

for ni=1:NImages
for nz=1:Nz
    temp_cmaps = squeeze(cmaps(nz,:,:,:));
    Images_raw(ni,nz,:,:) = squeeze(sum( squeeze(images_coils(nz,ni,:,:,:)) .* conj(temp_cmaps) ));    
end
end
clear ni nz temp_cmaps;

%% zero interpolation filling
if zero_params.onoff == 1
    Images = zeros(NImages, Nz, Nxy*zero_params.factor, Nxy*zero_params.factor); 
    for nz=1:Nz
        temp1 = squeeze(Images_raw(:,nz,:,:));
        if NImages==1
            temp2(1,:,:) = temp1;        
        else
            temp2 = temp1;
        end
        temp3 = mg_zero_filling(temp2, zero_params);
        if NImages==1
            temp4(1,:,:) = temp3;        
        else
            temp4 = temp3;
        end
        Images(:,nz,:,:) = temp4(:,:,:);
    end
else
    Images = Images_raw;
end
clear temp1 temp2 temp3 temp4 nz;

end
