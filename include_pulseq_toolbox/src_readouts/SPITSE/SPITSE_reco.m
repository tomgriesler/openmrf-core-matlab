function [Images, PULSEQ, study_info, cmaps, dcf2D, f0] = SPITSE_reco(study, cmaps, zero_params, mod_reco)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% ----- Input -----
% study:       enter path of study '.dat' or [] for dropdown select
% cmaps:       enter [] for calculating cmaps via openadapt or espirit
% zero_params: parameters for zero interpolation filling
% mod_reco:    1 -> openadapt()
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

%% import spiral rawdata
rawdata10D = twix_obj.image();
NRead      = size(rawdata10D,1);
NCoils     = size(rawdata10D,2);
NSeg       = PULSEQ.SPITSE.NEcho;
NImages    = size(rawdata10D,3)/NSeg;
f0         = study_info.f0;

%% sort segmented data
rawdata = zeros(NRead*NSeg, NCoils, NImages);
for j=1:NCoils
    temp1 = squeeze(rawdata10D(:,j,:));
    temp2 = temp1(:);
    if NImages>1
        for k=1:NImages
            temp3 = temp2( (k-1)*NRead*NSeg+1 : k*NRead*NSeg);
            rawdata(:,j,k) = temp3(:);
        end
    else
        rawdata(:,j,1) = temp2(:);
    end
    clear temp1 temp2;
end

rawdata(1:NRead/2,:,:) = [];

%% set nufft prams for reconstruction
PULSEQ.ktraj_reco(:,1:1:NRead/2) = []; % remove 1st prephaser
Nxy = PULSEQ.FOV.Nxy;
N   = [Nxy, Nxy];                      % [ ]   matrix size
fov = PULSEQ.FOV.fov_xy;               % [m]   fov size
kx  = PULSEQ.ktraj_reco(1,:,:);        % [1/m] x-trajecotry
ky  = PULSEQ.ktraj_reco(2,:,:);        % [1/m] y-trajecotry
kx  = kx(:);                           % join NRs
ky  = ky(:);                           % join NRs
kxy = [kx, ky];                        % Nadc x 2

%% calculate density compensation function: DCF
nufft_args = {N, [6 6], 2*N, N/2, 'table', 2^12, 'minmax:kb'};
nufft_st   = nufft_init( kxy / max(abs(kxy(:)))*pi, N, [6 6], N*2, N/2, 'minmax:kb');
G          = Gmri(kxy, true(N), 'fov', fov, 'basis', {'rect'}, 'nufft', nufft_args);
dcf        = abs(mri_density_comp( kxy, 'pipe', 'G', G.arg.Gnufft)) *2^2 *Nxy^2;
dcf2D      = mg_get_heatmap( kx, ky, dcf, 2 );

%% calculate nonuniform fast fourier transform: NUFFT
images_coils = zeros(NImages, NCoils, Nxy, Nxy);
parfor j = 1:NImages
    temp_k                = squeeze(rawdata(:,:,j));
    temp_reco             = permute( nufft_adj( temp_k .* repmat(dcf, [1, NCoils]), nufft_st ), [3, 1, 2]);
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

