function [rawdata_whitened, Psi, W] = mg_noise_prewhitening(rawdata, rawdata_noise, mode_w, flag_plot)

% INPUTS:
%   rawdata:       [Ncoils x ...] complex-valued raw MRI data
%   rawdata_noise: [Ncoils x Nnoise x Nadc] complex-valued noise-only data
%   mode_w:        'cholesky' or 'eigen'
%   flag_plot:     0 -> off; 1  -> on
%
% OUTPUT:
%   rawdata_whitened: [Ncoils x ...] prewhitened MRI data (same shape as input)
%   Psi:              [Ncoils x Ncoils] covariance matrix
%   W:                [Ncoils x Ncoils] whitening matrix

if nargin<3
    mode_w = [];
end
if nargin<4
    flag_plot = [];
end
if isempty(mode_w)
    mode_w = 'cholesky';
end
if isempty(flag_plot)
    flag_plot = 0;
end

% Reshape noise to 2D
Ncoils = size(rawdata_noise, 1);
rawdata_noise = reshape(rawdata_noise, Ncoils, []);  % [Ncoils x (Nnoise * Nadc)]

% Estimate noise covariance matrix
Psi = (rawdata_noise * rawdata_noise') / (size(rawdata_noise, 2)-1);  % [Ncoils x Ncoils]

% Compute whitening matrix: Cholesky factorization (option 1)
% http://hansenms.github.io/sunrise/sunrise2013/Hansen_ImageReconstruction_ParallelImaging.pdf
if strcmp(mode_w, 'cholesky')
    W = inv(chol(Psi, 'lower'));
end

% Compute whitening matrix: Eigenvalue decomposition (option 2)
if strcmp(mode_w, 'eigen')
    [U, Lambda] = eig(Psi);
    [lambda_sorted, idx] = sort(diag(Lambda), 'descend');
    U = U(:, idx);
    Lambda = diag(lambda_sorted);
    W = diag(1 ./ sqrt(diag(Lambda))) * U';  % [Ncoils x Ncoils]
end

% Normalize W to conserve noise energy
scale_factor = sqrt(trace(Psi) / trace(W * Psi * W'));
W = scale_factor * W;

% Reshape rawdata to 2D, apply whitening, reshape back
data_size  = size(rawdata);
num_total  = prod(data_size(2:end));               % total number of "time points" or spatial locations
rawdata_2D = reshape(rawdata, Ncoils, num_total);  % [Ncoils x num_voxels]
rawdata_whitened_2D = W * rawdata_2D;              % [Ncoils x num_voxels]
rawdata_whitened    = reshape(rawdata_whitened_2D, data_size);  % same shape as input

% vis covariance matrix and whitening matrix
if flag_plot==1
    figure()
    subplot(1,2,1)
    imagesc(abs(Psi)); axis image; colormap(hot); title('covariance matrix');
    subplot(1,2,2)
    imagesc(abs(W)); axis image; colormap(hot); title('whitening matrix');
end

end
