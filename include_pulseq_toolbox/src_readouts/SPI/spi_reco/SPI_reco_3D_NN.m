function [Images, PULSEQ, study_info] = SPI_reco_3D_NN(path_raw, path_backup, vendor, zero_params, ktraj_meas)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V2, 22.03.2026; unify for different vendors

% reconstruction for 3D datasets measured with, e.g. cones or seiffert spirals
% I don't need this, so I implemented a very simple
% nearest neighbour reconstruction without nufft!
% feel free to improve! :D

% case: no input arguments
if nargin==0
    path_raw    = [];
    path_backup = [];
    vendor      = [];
    zero_params = [];
    ktraj_meas  = [];
end

% defaults
if isempty(zero_params)
    zero_params.onoff  = 1;
    zero_params.radius = 6.0;
    zero_params.factor = 2.0;
end

% import rawdata, read study info and load pulseq workspace
[rawdata, rawdata_noise, PULSEQ, study_info] = pulseq_read_meas(path_raw, path_backup, vendor);

%% noise pre-whitening
if ~isempty(rawdata_noise)
    rawdata = mg_noise_prewhitening(rawdata, rawdata_noise, 'cholesky', 1);
end

%% svd coil compression
no_coils = 8;
if no_coils < size(rawdata,1)
    [~,~,v] = svd( reshape(rawdata,size(rawdata,1),[]).', 'econ');
    rawdata = permute( reshape(reshape(rawdata,size(rawdata,1),[]).'*v(:,1:no_coils), size(rawdata,2), size(rawdata,3), no_coils), [3,1,2]);
end

%% import measured k-space trajectories
if isempty(ktraj_meas)
    ktraj_reco = SPI_load_ktraj(PULSEQ);
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
end

%% Quick nearest-neighbour 3D gridding reconstruction: this was a single-prompt ChatGPT solution
N = PULSEQ.FOV.Nxy;

[Nc, NR, Nadc] = size(rawdata);
assert(all(size(ktraj_reco) == [3, NR, Nadc]), 'ktraj_reco must be 3 x NR x Nadc');

% --- reshape to sample list ---
S = NR * Nadc;
d = reshape(rawdata, [Nc, S]);                 % [Nc x S]
k = reshape(ktraj_reco, [3, S]);               % [3  x S]
kx = k(1,:); ky = k(2,:); kz = k(3,:);

% --- normalize k-space to grid indices 1..N (nearest neighbour) ---
% We estimate max radius across all axes (robust-ish quick scaling).
kmax = max([max(abs(kx)), max(abs(ky)), max(abs(kz))]);

% Map [-kmax, +kmax] -> [1, N]
ix = round( (kx./kmax * 0.5 + 0.5) * (N-1) + 1 );
iy = round( (ky./kmax * 0.5 + 0.5) * (N-1) + 1 );
iz = round( (kz./kmax * 0.5 + 0.5) * (N-1) + 1 );

% Keep only samples that fall inside the grid
mask = ix>=1 & ix<=N & iy>=1 & iy<=N & iz>=1 & iz<=N & isfinite(ix) & isfinite(iy) & isfinite(iz);
ix = ix(mask); iy = iy(mask); iz = iz(mask);
d  = d(:,mask);                                  % [Nc x S_in]

% Linear indices into 3D grid
lin = sub2ind([N N N], ix, iy, iz);              % [1 x S_in]

% Optional: average samples that hit the same Cartesian point (recommended)
cnt = accumarray(lin(:), 1, [N^3 1], @sum, 0);

% --- grid + IFFT coil-by-coil ---
img_coils = zeros(N, N, N, Nc, 'single');

for c = 1:Nc
    % Accumulate complex samples into Cartesian k-space grid (nearest neighbour)
    kgrid = accumarray(lin(:), double(d(c,:)).', [N^3 1], @sum, 0);
    
    % Average where multiple samples landed on same grid point
    kgrid(cnt > 0) = kgrid(cnt > 0) ./ cnt(cnt > 0);
    
    kgrid = reshape(kgrid, [N N N]);
    
    % 3D IFFT to image space (use shifts so k=0 is centered)
    img = ifftshift(ifftn(fftshift(kgrid)));
    
    img_coils(:,:,:,c) = single(img);
end

% --- simple coil combination ---
[Images, cmap] = openadapt(permute(img_coils,[4,1,2,3]));

%% zero interpolation filling
if zero_params.onoff == 1
    Images = mg_zero_filling(Images, zero_params);
end

end
