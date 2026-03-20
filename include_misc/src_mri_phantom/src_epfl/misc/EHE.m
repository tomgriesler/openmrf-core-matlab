%% EHE.m
%
% Function that performs the operation of the k-space measurement of a 2D
% MR experiment followed by its adjoint. Works for cartesian k-space
% trajectories without gridding.
%
% Inputs: the space domain original image, a structure containing the
% parameters (coils sensitivities, the k-space trajectory, convolution
% kernel generated e.g. by 'G_c.m').
%
% Output: a space domain (complex) image.
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group - EPFL, 2008-05-24

function data_out = EHE(data_in, param)

%% Initializations
res = size(data_in);
data_out = zeros(res);
Nb_coils = size(param.sensitivity,3);

%% Computation of the adjoint operator
% employ a precomputed convolution kernel
for ind_coil = 1:Nb_coils
    data = zeros(size(param.G));
    data(1:res(1),1:res(2)) = data_in.*param.sensitivity(:,:,ind_coil);
    data = fftn(data).*param.G;
    data = ifftn(data);
    data_out(:,:) = data_out(:,:) + conj(param.sensitivity(:,:,ind_coil)).*data(1:res(1),1:res(2));
end
data_out = data_out/prod(res);