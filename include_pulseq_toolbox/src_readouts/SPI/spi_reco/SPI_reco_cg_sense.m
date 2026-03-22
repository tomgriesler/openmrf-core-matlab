function [Images, PULSEQ, study_info, cmaps, dcf2D] = SPI_reco_cg_sense(path_raw, path_backup, vendor, cmaps, zero_params, n_iter, mod_reco, ktraj_meas)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V2, 22.03.2026; unify for different vendors

% ----- Input -----
% path_raw:    path of meas data or [] for select via uigetfile
% path_backup: path of pulseq workspace backup; not necessary for Siemens
% vendor:      vendor name
% cmaps:       enter [] for calculating cmaps via openadapt or espirit
% zero_params: parameters for zero interpolation filling
% mod_reco:    1 -> openadapt()
%              2 -> espirit()
% ktraj_meas:  [2 x Nid x Nadc] measured k-space trajectories															 

% case: no input arguments
if nargin==0
    path_raw    = [];
    path_backup = [];
    vendor      = [];
    cmaps       = [];
    zero_params = [];
    n_iter      = [];
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
if isempty(n_iter)
    n_iter = 25; % number of cg sense iterations
end

% import rawdata, read study info and load pulseq workspace
[rawdata, rawdata_noise, PULSEQ, study_info] = pulseq_read_meas(path_raw, path_backup, vendor);

%% read dimensions
NR = PULSEQ.SPI.NR;
[NCoils, NImages, NRead]  = size(rawdata);
NImages = NImages / NR;

%% noise pre-whitening
if ~isempty(rawdata_noise)
    rawdata = mg_noise_prewhitening(rawdata, rawdata_noise, 'cholesky', 1);
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
        rawdata   = rawdata(:,:, PULSEQ.SPI.adcNPad(1):PULSEQ.SPI.adcNPad(2));
        ktraj_reco = ktraj_reco(:,:, PULSEQ.SPI.adcNPad(1):PULSEQ.SPI.adcNPad(2));
    else
        rawdata   = rawdata(:,:,PULSEQ.SPI.adcNPad+1:end);
        ktraj_reco = ktraj_reco(:,:,PULSEQ.SPI.adcNPad+1:end);
    end
    NRead = size(rawdata, 3);
end

%% join NRs
rawdata2 = zeros(NCoils, NImages, NRead*NR);
for j = 1:NCoils
for k = 1:NImages
    temp = squeeze(rawdata(j, (k-1)*NR+1:k*NR ,:));
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

%% calculate coil sensitivity maps: cmaps
if isempty(cmaps)
    images_coils = zeros(NImages, NCoils, Nxy, Nxy);
    parfor j = 1:NImages
        temp_k                = squeeze(rawdata2(:,j,:));
        temp_reco             = permute( nufft_adj( temp_k.' .* repmat(dcf, [1, NCoils]), nufft_st ), [3, 1, 2]);
        images_coils(j,:,:,:) = temp_reco(:,:,:);
    end
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
    clear images_coils;
end

%% CG SENSE Recon

Images = zeros(NImages, Nxy, Nxy);
smaps  = permute(cmaps, [2,3,1]);
E      = @(x,tr) ismrm_encoding_non_cartesian_SENSE(x, smaps, nufft_st, dcf, tr);

parfor ni = 1:NImages

    temp_k   = squeeze(rawdata2(:,ni,:)).';
    reco_u   = nufft_adj( temp_k .* repmat(dcf *2^2*Nxy^2, [1 NCoils]), nufft_st );
    reco_u   = sum( conj(smaps) .* reco_u ,3) * sqrt(numel(dcf) / (2^4*Nxy^4));   
    a        = reco_u(:);
    b        = zeros(Nxy^2, 1);
    p        = a;
    r        = a;
    b_iter   = zeros(n_iter, Nxy, Nxy);
    d_iter   = zeros(n_iter, 1);
    
    for j = 1:n_iter
           d_iter(j,1) = (r' * r) / (a' * a);
           q  = E( E(p, 'notransp'), 'transp' );
           b  = b + (r' * r) / (p' * q) * p;      
           ri = r;       
           r  = ri - (ri' * ri) / (p'  * q)  * q;       
           p  = r  + (r'  * r)  / (ri' * ri) * p;
           b_iter(j,:,:) = reshape(b, [Nxy Nxy]);
    end

    Images(ni,:,:) = b_iter(end,:,:);

end

%% zero interpolation filling
if zero_params.onoff == 1
    Images = mg_zero_filling(Images, zero_params);
end

end
