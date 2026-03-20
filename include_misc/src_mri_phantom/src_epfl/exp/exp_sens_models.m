% Experiment 1: polynomial approximation of sensitivities
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 20-10-2009 (dd-mm-yyyy)
% Revised in early 2011

clear all;close all;
addpath('misc/');
make;

model_type = {'polynomial','sinusoidal'};
model_param_range = {[1:11],[1:9]};
model_plot = {'-rx','-bv'};
%Ntype = menu('Sensitivity model',model_type);
for Ntype = 1:length(model_type)
sens_model = model_type{Ntype};
param_range = model_param_range{Ntype};
plotopts = model_plot{Ntype};

%% Phantom
DefineSL;
res = 2^8*[1 1];
im = RasterizePhantom(SL,res);

%% Coils simulation
coil.Nb_coils = 24;
coil.res = res;
coil.type = 'biot';
coil.param.FS = 0.28; % FOV width
coil.param.D = 0.17; % Distance center->coil
coil.param.R = 0.05; % radius of the coil
coil.param.rand = 1;
coil = simulate_sensitivities(coil);

support = (im>1e-3);numel(find(support));
NbCoils = size(coil.sensitivity,3);
coil.sensitivity = coil.sensitivity/max(reshape(abs(coil.sensitivity.*repmat(support,[1,1,NbCoils])),1,numel(coil.sensitivity)));

%% Coils approximation

tic;
s=cell(1,NbCoils);
nrmse = zeros(length(param_range),NbCoils);
ser = zeros(length(param_range),NbCoils);
maxerror = zeros(length(param_range),NbCoils);
condi = zeros(length(param_range),1);
SER_mean = zeros(size(param_range));
NRMSE_mean = zeros(size(param_range));
maxerror_mean = zeros(size(param_range));
N_range = zeros(size(param_range));
for indD = 1:length(param_range)
    param = param_range(indD);
    for indCoil = 1:NbCoils
        [s{indCoil},nrmse(indD,indCoil),ser(indD,indCoil),maxerror(indD,indCoil),condi(indD)] = SensFitting(coil.sensitivity(:,:,indCoil),sens_model,param,support);
    end
    SER_mean(indD) = mean(ser(indD,:));
    NRMSE_mean(indD) = mean(nrmse(indD,:));
    maxerror_mean(indD) = mean(maxerror(indD,:));
    N_range(indD) = numel(s{1}.data);
    fprintf('%d coeffs . NRMSE : %g | relative max. error: %g\n',N_range(indD),NRMSE_mean(indD),maxerror_mean(indD));
end
toc
%%
%close all;
figure(1);hold on;plot(N_range,10*log10(SER_mean),plotopts);xlabel('number of parameters');title('Signal to Error ratio (dB)');hold off;
figure(2);hold on;semilogy(N_range,NRMSE_mean,plotopts);xlabel('number of parameters');set(gca,'Yscale','Log');title('NRMSE');hold off;
figure(3);hold on;semilogy(N_range,maxerror_mean,plotopts);xlabel('number of parameters');set(gca,'Yscale','Log');title('relative maximal absolute error');hold off;
figure(4);hold on;semilogy(N_range,condi,plotopts);xlabel('number of parameters');set(gca,'Yscale','Log');title('condition number of the fitting matrix');hold off;

%% Export results in LaTeX table
N_str = [];
NRMSE_str = [];
SER_str = [];
max_error_str = [];
dcl_str = '|c|';
for i = 1:length(N_range)
    dcl_str = [dcl_str, 'c|'];
    N_str = [N_str, sprintf(' %d',N_range(i)) ' &'];
    NRMSE_str = [NRMSE_str, sprintf(' %3.2e',NRMSE_mean(i)) ' &'];
    SER_str = [SER_str, sprintf(' %3.1f',10*log10(SER_mean(i))) ' &'];
    max_error_str = [max_error_str , sprintf(' %3.1e',maxerror_mean(i)) ' &'];
end
fexp1 = fopen(sprintf('../table_exp1_%s.tex',sens_model),'w+');
fprintf(fexp1,'\\begin{center}\n\\begin{tabular}{%s}\n\\hline Number of coefficients & %s \\\\\n\\hline\\hline NRMSE & %s \\\\\n\\hline SER (dB) & %s \\\\\n\\hline max.\\ error & %s \\\\\n\\hline\n\\end{tabular}\n\\end{center}',dcl_str,N_str(1:end-1),NRMSE_str(1:end-1),SER_str(1:end-1),max_error_str(1:end-1));
fclose(fexp1);
end
%%
figure(1);legend(model_type);
figure(2);legend(model_type);
figure(3);legend(model_type);
figure(4);legend(model_type);
