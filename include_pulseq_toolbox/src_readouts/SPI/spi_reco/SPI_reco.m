function [Images, PULSEQ, study_info, cmaps, dcf2D, f0, rawdata_noise] = SPI_reco(study, cmaps, zero_params, mod_reco, ktraj_meas)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% ----- Input -----
% study:       enter path of study '.dat' or [] for dropdown select
% cmaps:       enter [] for calculating cmaps via openadapt or espirit
% zero_params: parameters for zero interpolation filling
% mod_reco:    1 -> openadapt()
%              2 -> espirit()
% ktraj_meas:  [2 x Nid x Nadc] measured k-space trajectories

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
NImages       = NImages / NR;
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

%% join NRs
rawdata2 = zeros(NCoils, NImages, NRead*NR);
for j = 1:NCoils
for k = 1:NImages
    temp = squeeze(rawdata1(j, (k-1)*NR+1:k*NR ,:));
    temp = temp(:);
    rawdata2(j, k, :) = temp(:);
end
end
clear j k temp;

%% set nufft prams for reconstruction
Nxy = PULSEQ.FOV.Nxy;
N   = [Nxy, Nxy];          % [ ]   matrix size
fov = PULSEQ.FOV.fov_xy;   % [m]   fov size
kx  = ktraj_reco(1,:,:);   % [1/m] x-trajecotry
ky  = ktraj_reco(2,:,:);   % [1/m] y-trajecotry
kx  = kx(:);               % join NRs
ky  = ky(:);               % join NRs
kxy = [kx, ky];            % Nadc x 2

%% calculate density compensation function: DCF
nufft_args = {N, [6 6], 2*N, N/2, 'table', 2^12, 'minmax:kb'};
nufft_st   = nufft_init( kxy / max(abs(kxy(:)))*pi, N, [6 6], N*2, N/2, 'minmax:kb');
G          = Gmri(kxy, true(N), 'fov', fov, 'basis', {'rect'}, 'nufft', nufft_args);
dcf        = abs(mri_density_comp( kxy, 'pipe', 'G', G.arg.Gnufft)) *2^2 *Nxy^2;
dcf2D      = mg_get_heatmap( kx, ky, dcf, 2 );

%% calculate nonuniform fast fourier transform: NUFFT
images_coils = zeros(NImages, NCoils, Nxy, Nxy);
parfor j = 1:NImages
    temp_k                = squeeze(rawdata2(:,j,:));
    temp_reco             = permute( nufft_adj( temp_k.' .* repmat(dcf, [1, NCoils]), nufft_st ), [3, 1, 2]);
    images_coils(j,:,:,:) = temp_reco(:,:,:);
end

%% calculate coil sensitivity maps: cmaps
if isempty(cmaps)
    if mod_reco==1
        if NImages>1
            [~, cmaps] = openadapt(squeeze(mean(images_coils)));
        else
            [~, cmaps] = openadapt(squeeze(images_coils));
        end
    end
    if mod_reco==2
        if NImages>1
            cmaps = mg_espirit_cmaps(squeeze(mean(images_coils)), 0.02, 0.95, 24, [6,6]);
        else
            cmaps = mg_espirit_cmaps(squeeze(images_coils), 0.02, 0.95, 24, [6,6]);
        end
    end
end

%% calculate coil combined images
Images = zeros(NImages, Nxy, Nxy);
for j=1:NImages 
    Images(j,:,:) = squeeze(sum( squeeze(images_coils(j,:,:,:)) .* conj(cmaps) ));    
end

%% zero interpolation filling
if zero_params.onoff == 1
    Images = mg_zero_filling(Images, zero_params);
end

end
