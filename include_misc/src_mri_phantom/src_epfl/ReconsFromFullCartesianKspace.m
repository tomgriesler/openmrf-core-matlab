%% ReconsFromFullCartesianKspace.m
%
% Reconstruc in space MR data that where acquired on a full Cartesian
% kspace grid.
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 31-10-2009 (dd-mm-yyyy)
% revised in Feb. 2011

function im = ReconsFromFullCartesianKspace(m,k,FOV)

k = diag(FOV)*k;
res = ceil(max(k,[],2)-min(k,[],2)+1)';

im = zeros(res);
idx1 = mod(round(-k(1,:)),res(1))+1;
idx2 = mod(round(-k(2,:)),res(2))+1;
ind = sub2ind(res,idx1,idx2);
im(ind) = m;%imagesc(log10(abs(im)));pause
im = fftshift(ifftn(im))*prod(res)/prod(FOV);