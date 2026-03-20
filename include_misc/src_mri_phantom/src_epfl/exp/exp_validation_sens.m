close all;clc;

res = 128;

%% PHANTOM

DefineBrain;
DefineSL;
leg = {'brain', 'Sheep-Logan'};

%% SENSITIVITY

%% Coils simulation
coil.Nb_coils = 1;
coil.res = res;
coil.type = 'biot';
coil.param.FS = 0.28; % FOV width
coil.param.D = 0.17; % Distance center->coil
coil.param.R = 0.05; % radius of the coil
coil.param.rand = 0;
coil = simulate_sensitivities(coil);

NbCoils = size(coil.sensitivity,3);

sens.model = 'sinusoidal';
%sens.model = 'polynomial';
sens.param = 7;

im = RasterizePhantom(Brain,res,[1],0);
support = (im>1e-3);numel(find(support));
sensitivity = coil.sensitivity/max(reshape(abs(coil.sensitivity.*repmat(support,[1,1,NbCoils])),1,numel(coil.sensitivity)));
sens_brain = SensFitting(sensitivity(:,:,1),sens.model,sens.param,support);

im = RasterizePhantom(SL,res,[1],0);
support = (im>1e-3);numel(find(support));
sensitivity = coil.sensitivity/max(reshape(abs(coil.sensitivity.*repmat(support,[1,1,NbCoils])),1,numel(coil.sensitivity)));
sens_sl = SensFitting(sensitivity(:,:,1),sens.model,sens.param,support);

% if 0%HOMOGENEOUS
%     sens.data = 1;
%     sens.param = 1;
% else
%     sens.data = [0,-.2,0,.2,.5,.2,0,-.2,0]';
%     sens.param = 3;
% end

%%
res_range = [128,176,256,352,512,704];%,1024,1408,2048,2816

%% K-SPACE
k = GenerateFullCart2DKspace(128*[1,1]);

%% MR DATA
tic;m_brain = MRData(Brain,sens_brain,k);t_analytical(1)=toc;
tic;m_sl = MRData(SL,sens_sl,k);t_analytical(2)=toc;
%%
rec_brain = ReconsFromFullCartesianKspace(m_brain,k,Brain.FOV);
rec_sl = ReconsFromFullCartesianKspace(m_sl,k,Brain.FOV);
for i = 1:length(res_range)
    res= res_range(i)*[1,1];
    
    tic;
    im_raster = RasterizePhantom(Brain,res,sens_brain,0);
    m_num = MRDataNumerical(im_raster,k,Brain.FOV);
    t_num(i,1)=toc;
    rec_num = ReconsFromFullCartesianKspace(m_num,k,Brain.FOV);
    snr(i,1) = norm(rec_num(:),2)/norm(rec_num(:)-rec_brain(:),2);
    maxerr_space(i,1) = max(abs(rec_num(:)-rec_brain(:)))/max(abs(rec_num(:)));
    maxerr_kspace(i,1) = max(abs(m_num(:)-m_brain(:)))/max(abs(m_brain(:)));
    tic;
    im_raster = RasterizePhantom(SL,res,sens_sl,0);
    m_num = MRDataNumerical(im_raster,k,SL.FOV);
    t_num(i,2)=toc;
    rec_num = ReconsFromFullCartesianKspace(m_num,k,SL.FOV);
    snr(i,2) = norm(rec_num(:),2)/norm(rec_num(:)-rec_sl(:),2);
    maxerr_space(i,2) = max(abs(rec_num(:)-rec_sl(:)))/max(abs(rec_num(:)));
    maxerr_kspace(i,2) = max(abs(m_num(:)-m_sl(:)))/max(abs(m_sl(:)));
    
    %
    %imax = max([abs(rec_num(:));abs(rec_(:))]);
    
    figure(1);semilogx(res_range(1:i),20*log10(snr(1:i,1)),'b.-',res_range(1:i),20*log10(snr(1:i,2)),'r.-');ylabel('SNR');xlabel('sampling density');legend(leg)
    figure(5);loglog(res_range(1:i),1./snr(1:i,1),'b.-',res_range(1:i),1./snr(1:i,2),'r.-');ylabel('NRMSE');xlabel('sampling density');legend(leg)
    figure(2);loglog(res_range(1:i),maxerr_space(1:i,1),'b.-',res_range(1:i),maxerr_space(1:i,2),'r.-');ylabel('max error in space');xlabel('sampling density');legend(leg)
    figure(3);loglog(res_range(1:i),maxerr_kspace(1:i,1),'b.-',res_range(1:i),maxerr_kspace(1:i,2),'r.-');ylabel('max error in kspace');xlabel('sampling density');legend(leg)
    figure(4);loglog(res_range(1:i),t_num(1:i,1)/t_analytical(1),'b.-',res_range(1:i),t_num(1:i,2)/t_analytical(2),'r.-');ylabel('t_{num}/t_{analytical}');xlabel('sampling density');legend(leg)
end
%save 'exp_validation_sens.mat' leg res_range t_analytical t_num maxerr_kspace maxerr_space snr;
%%
%figure('Name','Rasterized');imagesc(abs(rec_num),[0,imax]);axis image;colorbar;%colormap gray;
%figure('Name','Reconstructed');imagesc(abs(rec),[0,imax]);axis image;colorbar;%colormap gray;
%figure('Name',sprintf('Diff with rasterized. SNR = %f',20*log10(norm(rec_num(:),2)/norm(rec_num(:)-rec(:),2))));imagesc(abs(rec_num-rec),[0,imax]);axis image;colorbar;%colormap gray;
% clear all;
% load exp_validation_sens.mat;

N_str = [];
NRMSE_brain_str = [];
max_errork_brain_str = [];
max_error_brain_str = [];
t_brain_str = [];
NRMSE_sl_str = [];
max_errork_sl_str = [];
max_error_sl_str = [];
t_sl_str = [];
dcl_str = '|c|c|';
for i = 1:length(res_range)
    dcl_str = [dcl_str, 'c|'];
    N_str = [N_str, sprintf(' %d\t\t',res_range(i)) '&'];
    NRMSE_brain_str = [NRMSE_brain_str, sprintf(' %3.2e\t',1/snr(i,1)) '&'];
    max_errork_brain_str = [max_errork_brain_str, sprintf(' %3.1e\t',maxerr_kspace(i,1)) '&'];
    max_error_brain_str = [max_error_brain_str , sprintf(' %3.1e\t',maxerr_space(i,1)) '&'];
    t_brain_str = [t_brain_str, sprintf(' %.4f\t',t_num(i,1)/t_analytical(1)) '&'];
    NRMSE_sl_str = [NRMSE_sl_str, sprintf(' %3.2e\t',1/snr(i,2)) '&'];
    max_errork_sl_str = [max_errork_sl_str, sprintf(' %3.1e\t',maxerr_kspace(i,2)) '&'];
    max_error_sl_str = [max_error_sl_str , sprintf(' %3.1e\t',maxerr_space(i,2)) '&'];
    t_sl_str = [t_sl_str, sprintf(' %.4f\t',t_num(i,2)/t_analytical(2)) '&'];
end
sprintf('\\begin{center}\n\\begin{tabular}{%s}\n\\hline\\hline\n \t& Sampling density\t& %s \\\\\n\\hline \t& NRMSE\t\t\t& %s \\\\\n Brain \t& max. err. k-space\t& %s \\\\\n \t& max. err. pixels\t& %s \\\\\n \t& rel. time\t\t& %s \\\\\n\\hline\n \t& NRMSE\t\t\t& %s \\\\\n SL \t& max. err. k-space\t& %s \\\\\n \t& max. err. pixels\t& %s \\\\\n \t& rel. time\t\t& %s \\\\\n\\hline\\hline\n\\end{tabular}\n\\end{center}',dcl_str(2:end-1),N_str(1:end-1),NRMSE_brain_str(1:end-1),max_errork_brain_str(1:end-1),max_error_brain_str(1:end-1),t_brain_str(1:end-1),NRMSE_sl_str(1:end-1),max_errork_sl_str(1:end-1),max_error_sl_str(1:end-1),t_sl_str(1:end-1))
