%% GenerateFullCart2DKspace.m
%
% Function that generates a 2D full Cartesian grid.
%
% INPUT:    * res : (2,1) vector or scalar that defines resolution
%
% OUTPUT:   * k : (2,N) matrix of kspace trajectory
% 
% Note: N = nb of samples = prod(res)
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 30-10-2009 (dd-mm-yyyy)

function k = GenerateFullCart2DKspace(res,FOV)

if numel(res)==1
    res = res*[1,1];
end
if nargin<2
    FOV = [1,1];
end

[kX,kY] = GenerateFullCart2DGrid(res);

k = diag(1./FOV)*[kX(:)';kY(:)'];
