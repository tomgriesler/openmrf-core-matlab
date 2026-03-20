function [cmaps, mask, im] = mg_espirit_cmaps(images_coils, eigThresh_1, eigThresh_2, ncalib, ksize)

% mgram: code was copied from demo_ESPIRiT_maps.m

% ----- input -----
% images_coils: fully sampled images -> Ncoils x Ny x Nx
% eigThresh_1:  threshold for picking singular vercors of the calibration matrix
% eigThresh_2:  threshold of eigen vector decomposition in image space
% ncalib:       use 24 calibration lines to compute compression
% ksize:        kernel size
% ----- output -----
% cmaps: coil sensitivity maps -> Ncoils x Ny x Nx
% im:    combined image -> Ny x Nx


if isempty(eigThresh_1)
    eigThresh_1 = 0.02;
end
if isempty(eigThresh_2)
    eigThresh_2 = 0.95;
end
if isempty(ncalib)
    ncalib = 24;
end
if isempty(ksize)
    ksize = [6, 6];
end

[Nc, sy, sx] = size(images_coils);
kspace_coils = image2kspace(images_coils);
kspace_coils = permute(kspace_coils, [2,3,1]);

% crop a calibration area
calib = crop(kspace_coils, [ncalib, ncalib, Nc]);

% compute Calibration matrix
% perform 1st SVD and convert singular vectors
% into k-space kernels
[k, S] = dat2Kernel(calib,ksize);
idx    = max(find(S >= S(1)*eigThresh_1));
    
% crop kernels and compute eigen-value decomposition in image space
[M, W] = kernelEig(k(:,:,:,1:idx),[sy,sx]);

% crop sensitivity maps
mask  = W(:,:,end)>eigThresh_2;
cmaps = M(:,:,:,end);% .* repmat(mask, [1,1,Nc]); % mgram: get cmaps without mask
cmaps = permute(cmaps, [3,1,2]);

% ombine coil signals
im = squeeze(sum(images_coils .* conj(cmaps)));

end

