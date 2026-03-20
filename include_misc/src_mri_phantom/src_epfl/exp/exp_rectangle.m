close all;clc;
res = 256*[1,1];

%% PHANTOM
DefineRect;

%% SENSITIVITY
HOMOGENEOUS = true;
sens.model = 'polynomial';
if HOMOGENEOUS
    sens.data = 1;
    sens.param = 0;
else
    sens.data = [1,.2,.2,-2,.1,5]';
    sens.param = 2; 
end

% sens.model = 'sinusoidal';
% if HOMOGENEOUS
%     sens.data = 1;
%     sens.param = 1;
% else
%     sens.data = [0,-.2,0,.2,.5,.2,0,-.2,0]';
%     sens.param = 3; 
% end

%% K-SPACE
k = GenerateFullCart2DKspace(res);

%% MR DATA

tic;m = MRData(rect,sens,k);toc
rec = ReconsFromFullCartesianKspace(m,k,rect.FOV);

if HOMOGENEOUS
    m_theoretical = prod(sinc(diag(width)*R*k),1)*prod(width);
    fprintf('K-space:\n\t-> max error = %g,\n\t-> mean error = %g\n',max(abs(m_theoretical(:)-m(:)))/max(abs(m(:))),norm(m_theoretical(:)-m(:),2)/norm(m_theoretical(:),2));
    rec_theoretical = ReconsFromFullCartesianKspace(m_theoretical,k,rect.FOV);
    fprintf('Space:\n\t-> max error = %g,\n\t-> mean error = %g\n',max(abs(rec_theoretical(:)-rec(:)))/max(abs(rec_theoretical(:))),norm(rec_theoretical(:)-rec(:),2)/norm(rec_theoretical(:),2));
    %figure('Name',sprintf('Diff with theoretical. SNR = %f',20*log10(norm(rec_theoretical(:),2)/norm(rec_theoretical(:)-rec(:),2))));imagesc(abs(rec_theoretical-rec),[0,imax]);axis image;%colorbar;
end

res_range = [256,352,400,512,704,800,1024,1408,1600,2048];
for i=1:length(res_range)
    res= res_range(i)*[1,1];
    im_raster = RasterizePhantom(rect,res,sens,0);
    m_num = MRDataNumerical(im_raster,k,rect.FOV);
    rec_num = ReconsFromFullCartesianKspace(m_num,k,rect.FOV);
    snr(i) = norm(rec_num(:),2)/norm(rec_num(:)-rec(:),2);
    maxerr_space(i) = max(abs(rec_num(:)-rec(:)))/max(abs(rec_num(:)));
    maxerr_kspace(i) = max(abs(m_num(:)-m(:)))/max(abs(m(:)));
    
%imax = max([abs(im_raster(:));abs(rec(:))]);
%figure('Name','Rasterized');imagesc(abs(im_raster),[0,imax]);axis image;colorbar;%colormap gray;
%figure('Name','Reconstructed form analytical');imagesc(abs(rec),[0,imax]);axis image;colorbar;%colormap gray;
%figure('Name',sprintf('Diff with rasterized. SNR = %f',20*log10(norm(im_raster(:),2)/norm(im_raster(:)-rec(:),2))));imagesc(abs(im_raster-rec),[0,imax]);axis image;colorbar;%colormap gray;
end
figure(1);semilogx(res_range,20*log10(snr),'b.-');ylabel('SNR');xlabel('sampling density');
figure(5);loglog(res_range,1./snr,'b.-');ylabel('NRMSE');xlabel('sampling density');
figure(2);loglog(res_range,maxerr_space,'b.-');ylabel('max error in space');xlabel('sampling density');
figure(3);loglog(res_range,maxerr_kspace,'b.-');ylabel('max error in kspace');xlabel('sampling density');


%%
N_str = [];
NRMSE_str = [];
max_errork_str = [];
max_error_str = [];
dcl_str = '|l|';
for i = 1:length(res_range)
    dcl_str = [dcl_str, 'c|'];
    N_str = [N_str, sprintf(' %d',res_range(i)) ' &'];
    NRMSE_str = [NRMSE_str, sprintf(' %3.2e',1/snr(i)) ' &'];
    max_errork_str = [max_errork_str, sprintf(' %3.1e',maxerr_kspace(i)) ' &'];
    max_error_str = [max_error_str , sprintf(' %3.1e',maxerr_space(i)) ' &'];
end
sprintf('\\begin{center}\n\\begin{tabular}{%s}\n\\hline\\hline Sampling density & %s \\\\\n\\hline NRMSE & %s \\\\\n max. err. k-space & %s \\\\\n max. err. pixels & %s \\\\\n\\hline\\hline\n\\end{tabular}\n\\end{center}',dcl_str(2:end-1),N_str(1:end-1),NRMSE_str(1:end-1),max_errork_str(1:end-1),max_error_str(1:end-1))
