%% MRDataNumerical.m
%
% Function that generates numerically some kspace data given rasterized
% spatial data.
%
% To reduce aliasing artifacts, the resolution of the spatial data should
% be as large as possible.
%
% This code works for non-Cartesian kspace samples thanks to the NUFFT
% routines of Fessler and Sutton. The 'irt' toolbox must be installed.
% SEE: http://www.eecs.umich.edu/~fessler/code/
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 31-10-2009 (dd-mm-yyyy)
% revised in Feb. 2011

function m = MRDataNumerical(pix_phant,k,FOV)

k = diag(FOV)*k;
ref = 1; % refinement
tol = 1e-8;
rs = size(pix_phant);
if norm(round(ref*k(:))-ref*k(:),2)<tol % in case of Cartesian sampling
    %disp('fft');
    X = zeros(ref*rs);
    X(1:rs(1),1:rs(2)) = pix_phant;
    X = circshift(X,-floor(rs/2));%imagesc(abs(X));pause
    X = fftn(X);
    ind = sub2ind(ref*rs,mod(round(-ref*k(1,:)),ref*rs(1))+1,mod(round(-ref*k(2,:)),ref*rs(2))+1);
    m = X(ind);
else % for non-Cartesian sampling
    %disp('nufft');
    st = nufft_init(-2*pi*[k(1,:)'/rs(1),k(2,:)'/rs(2)], rs, [10,10], 2*rs);
    %m = nufft(circshift(pix_phant,-floor(rs/2)), st);
    m = nufft(pix_phant, st);
    m = (m.').*exp(-1i*pi*(k(1,:)+k(2,:)));
end
m = prod(FOV)/prod(rs)*m;
%Note: circshift is used to count the origin as the center of the image
% alternatively m can be multiplied by .*exp(-pi*1j*(k(1,:)+k(2,:)));
