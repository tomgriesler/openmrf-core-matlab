close all;clc;

res = 64;

%% PHANTOM

DefineSL;

%% Coils simulation
coil.Nb_coils = 1;
coil.res = res;
coil.type = 'homogeneous';
coil = simulate_sensitivities(coil);

NbCoils = size(coil.sensitivity,3);

sens.model = 'sinusoidal';
sens.data = 1;
sens.param = 1;

%% K-SPACE
k = GenerateFullCart2DKspace(res*[1,1]);

%% MR DATA
tic;m_analytical = MRData(SL,sens,k);t_analytical=toc;
rec_analytical = ReconsFromFullCartesianKspace(m_analytical,k,SL.FOV);
im = RasterizePhantom(SL,res,sens,0);
tic;m_raster = MRDataNumerical(im,k,SL.FOV);t_num=toc;
rec_raster = ReconsFromFullCartesianKspace(m_raster,k,SL.FOV);

%%
close all;

imax = max(max(real(rec_analytical(:))),max(real(rec_raster(:))));
imin = min(min(real(rec_analytical(:))),min(real(rec_raster(:))));
dyn = imax-imin;

x = res/2+1;

figure('Name','ANALYTICAL');imagesc(real(rec_analytical-imin)/dyn);axis image;colormap gray;
figure('Name','RASTERIZED');imagesc(real(rec_raster-imin)/dyn);axis image;colormap gray;
figure('Name','DIFFERENCE');imagesc(abs(rec_analytical-rec_raster)/dyn);axis image;colormap gray;

figure('Name','Profile');plot(1:res,real(rec_raster(x,:)),'g-*','LineWidth',1)
hold on;plot(1:res,real(rec_analytical(x,:)),'r+-','LineWidth',1);%,1:res,imag(rec_analytical(res/2,:)),'gx');
%legend('Rasterized','Real part of analytical','Imaginary part of analytical');
legend('Rasterized','Analytical');

print_analyt = repmat(uint16(2^16*real(rec_analytical-imin)/dyn),[1,1,3]);
print_analyt(x,:,1) = 2^16-1;
print_analyt(x,:,2) = 0;
print_analyt(x,:,3) = 0;
print_raster = repmat(uint16(2^16*real(rec_raster-imin)/dyn),[1,1,3]);
print_raster(x,:,1) = 0;
print_raster(x,:,2) = 2^16-1;
print_raster(x,:,3) = 0;

imwrite(print_raster,'rasterized.png');
imwrite(print_analyt,'analytical.png');