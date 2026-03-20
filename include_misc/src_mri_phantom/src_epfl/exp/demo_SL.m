%% Demonstration of analytical computations with the 2D SL brain phantom
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 20-10-2009 (dd-mm-yyyy)

clear all;close all;clc
addpath('misc/');
make;
gamma = 1;

model_type = {'polynomial','exponential'};
model_param = {8,6};
Ntype = menu('Sensitivity model',model_type);
sens_model = model_type{Ntype};
param = model_param{Ntype};

%% Phantom
DefineSL;
FOV = SL.FOV;
res = 2^9*[1,1];
h1 = figure('Name','Rasterized phantom','NumberTitle','off');
phant = RasterizePhantom(SL,res,1);imagesc(phant);axis image;colormap gray;

%% K-space trajectory
res_traj = 2^8*[1 1];

% generate k-space trajectory (full Cartesian in that case)
k = GenerateFullCart2DKspace(res_traj,FOV);

%% Coils

coil.Nb_coils = 1;
coil.res = res;
coil.type = 'biot';
coil.param.FS = 0.28; % FOV width
coil.param.D = 0.25; % Distance center->coil
coil.param.R = 0.1; % radius of the coil
coil = simulate_sensitivities(coil);
%%

support = (phant>1e-3);
h5 = figure('Name','Numerical computations','NumberTitle','off');
s = cell(1,size(coil.sensitivity,3));
sens_num = cell(1,size(coil.sensitivity,3));
residue = cell(1,size(coil.sensitivity,3));
tic;
for i = 1:size(coil.sensitivity,3)
    sens = coil.sensitivity(:,:,i);
    sens = sens/max(sens(support));
    [s{i},nrmse,ser,maxerrorsens] = SensFitting(sens,sens_model,param,support);
    [im,sens_num{i}] = RasterizePhantom(SL,res,s{i});
    residue{i} = sens_num{i}(support)-sens(support);
    subplot(1,size(coil.sensitivity,3),i);imagesc(abs(im).^gamma);colormap gray;axis image;%colorbar;
    fprintf('coil %d . Signal to error ratio (dB): %g\n max error: %g\n',i,10*log10((sens(support)'*sens(support))./(residue{i}(:)'*residue{i}(:))),max(abs(residue{i}(:)))/max(abs(sens(support))));
end
toc

%% Computations
m = cell(1,length(s));
h2 = figure('Name','Analytical Computations','NumberTitle','off');
h3 = figure('Name','Differences with numerical','NumberTitle','off');
rec = cell(1,length(s));
difference = cell(1,length(s));
m_num = cell(1,length(s));
rec_num = cell(1,length(s));
tic;
for i = 1:length(s)
    m{i} = MRData(SL,s{i},k);
    rec{i} = ReconsFromFullCartesianKspace(m{i},k,FOV);    
    m_num{i} = MRDataNumerical(sens_num{i}.*phant,k,FOV);
    rec_num{i} = ReconsFromFullCartesianKspace(m_num{i},k,FOV);    
    difference{i} = rec{i}-rec_num{i};    
    figure(h2);subplot(1,length(s),i);imagesc((abs(rec{i})).^gamma);colormap gray;axis image;%colorbar;title('computations in space');
    figure(h3);subplot(1,length(s),i);imagesc((abs(difference{i})).^gamma);colormap gray;axis image;%colorbar;title('computations in space');
end
toc

%% Save images
% fshep = fopen(sprintf('../fig_shepp_%s.tex',sens_model),'w+');
% fprintf(fshep,'\\begin{center}\n');
% for i = 1:length(param)
%     imwrite(uint8(round(255*(abs(rec{i}).^gamma)/max((abs(rec{i}(:)).^gamma)))),['../fig/shepp' num2str(i) '.png']);
%     fprintf(fshep,'\\includegraphics[width=%0.3f\\textwidth]{shepp%d.png}\n',0.95/length(param),i);
%     imwrite(uint8(round(255*(abs(difference{i}).^gamma)/max(abs(rec{i}(:)).^gamma))),['../fig/shepp_diff_' num2str(i) '.png']);
% end
% fprintf(fshep,'\\\n');
% for i = 1:length(param)
%     fprintf(fshep,'\\includegraphics[width=%0.3f\\textwidth]{shepp_diff_%d.png}\n',0.5*0.95/length(param),i);
% end
% for i = 1:length(s)
%     imwrite(uint8(round(255*(abs(rec{i})/max(abs(rec{i}(:)))).^gamma)),['../fig/shepp' num2str(i) '.png']);
%     fprintf(fshep,'\\includegraphics[width=%0.3f\\textwidth]{shepp%d.png}\\hspace{0.2cm}\\includegraphics[width=%0.3f\\textwidth]{shepp_diff_%d.png}\\\\\n\\vspace{0.2cm}\n',0.90*0.25,i,0.90*0.25,i);
%     imwrite(uint8(round(255*(abs(difference{i}/max(abs(difference{i}(:))))).^gamma)),['../fig/shepp_diff_' num2str(i) '.png']);
% end
% 
% fprintf(fshep,'\\end{center}\n');
% fclose(fshep);