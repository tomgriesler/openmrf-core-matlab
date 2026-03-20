%% Experiment 3: Computation errors with rasterization-based simulation
%               + interpolation errors
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 20-10-2009 (dd-mm-yyyy)

clear all;close all;
addpath('misc/');
%make;

model_type = {'polynomial','sinusoidal'};
model_param_range = {[8],[6]};
Ntype = menu('Sensitivity model',model_type);
sens_model = model_type{Ntype};
param_range = model_param_range{Ntype};

%% Phantom
DefineSL;
FOV = SL.FOV;

%% Coils

res_ref = 128*[1 1];

coil.Nb_coils = 3;
coil.res = res_ref;
coil.type = 'biot';%'homogeneous';%
coil.param.FS = 0.28; % FOV width
coil.param.D = 0.20; % Distance center->coil
coil.param.R = 0.05; % radius of the coil
%coil.param.rand = 1;
coil = simulate_sensitivities(coil);

%% Polynomial approximation of coils

l_range = 1:3;
color = {'b','c','r','m','g','y','k'};
phant = RasterizePhantom(SL,res_ref);
support = (phant>1e-3);
t_ref = zeros(1,length(param_range));
t_discrete = cell(1,length(param_range));
NRMSE = cell(1,length(param_range));
SER = cell(1,length(param_range));
maxerror = cell(1,length(param_range));
leg = cell(1,length(param_range));

for c = 1:size(coil.sensitivity,3)
    sens = coil.sensitivity(:,:,c);%clear coil;
    for b = 1:length(param_range)
        
        param = param_range(b);        
        [s,mse,ser,maxerrorsens] = SensFitting(sens,sens_model,param,support);
        
        %% K-space trajectory
        res_traj = res_ref;
        % Get k-space trajectory
        load spiral_traj;
        k = [k(1,:)*res_traj(1)/FOV(1);k(2,:)*res_traj(2)/FOV(2)]*0.9;
        
        %% perform analytical computations
        t0 = clock();
        m = MRData(SL,s,k);
        t_ref(b) = etime(clock(),t0);
        
        
        %% numerical computations
        t_discrete{b} = zeros(1,length(l_range));
        tic;
        for l = 1:length(l_range)
            res = 2^l_range(l)*res_ref;
            
            t0 = clock();
            % generate the rasterized phantom weighted by the sensitivity
            rasterized_weighted_phant = RasterizePhantom(SL,res,s);
            
            % discretized simulations
            m_num = MRDataNumerical(rasterized_weighted_phant,k,FOV);
            t_discrete{b}(l) = etime(clock(),t0);
            % error
            Error = m-m_num;
            MSE = Error(:)'*Error(:)/numel(Error);
            NRMSE{b}(l) = sqrt(MSE)/(max(abs(m(:))-min(abs(m(:)))));
            SER{b}(l) = (m(:)'*m(:))/MSE/numel(m);
            maxerror{b}(l) = max(abs(m(:)-m_num(:)))/(max(abs(m(:))-min(abs(m(:)))));
            fprintf('size %dx%d . NRMSE : %3.2e | SER (dB) : %3.2f | relative max. error: %3.2e \n',res(1),res(2),NRMSE{b}(l),10*log10(SER{b}(l)),maxerror{b}(l));
        end
        toc
        
        counter = 1;
        tick = cell(1,length(l_range));
        for i=l_range(1:end)
            tick{counter} = [num2str(2^l_range(i)*res_ref(1)) 'x' num2str(2^l_range(i)*res_ref(2))];
            counter = counter+1;
        end
        leg{b} = ['sensitivity degree: ' num2str(param_range(b))];
        figure(12);
        hold on;
        plot(l_range,10*log10(SER{b}),['*-' color{mod(c-1,length(color))+1}]);xlabel('resolution of rasterization');ylabel('Signal to error ratio (dB)');
        %legend(leg);hold off;
        
        figure(13);
        hold on;
        semilogy(l_range,t_discrete{b}/t_ref(b),['*-' color{mod(c-1,length(color))+1}]);xlabel('resolution of rasterization');ylabel('speedup');
        %legend(leg);hold off;
        
        figure(14);
        hold on;
        semilogy(l_range,maxerror{b},['*-' color{mod(c-1,length(color))+1}]);xlabel('resolution of rasterization');ylabel('max. error');
    end
end
figure(12);hold on;
plot([min(l_range),max(l_range)],30*[1,1],'k-','LineWidth',2);hold off;
set(gca,'XTick',l_range(1:end));
set(gca,'XTickLabel',tick);
figure(13);hold on;
plot([min(l_range),max(l_range)],1*[1,1],'k-','LineWidth',2);hold off;
set(gca,'XTick',l_range(1:end));
set(gca,'XTickLabel',tick);
figure(14);
set(gca,'XTick',l_range(1:end));
set(gca,'XTickLabel',tick);

%% Latex table

%upsmpl_str = [];
nb_pix_str = [];
NRMSE_str = [];
SER_str = [];
maxerror_str = [];
speedup_str = [];
time_str = [];
dcl_str = '|c|';
for i = 1:length(l_range)
    dcl_str = [dcl_str, 'c|'];
    %upsmpl_str = [upsmpl_str, sprintf(' %d',l_range(i)-log2(min(res_traj))-1) ' &'];
    nb_pix_str = [nb_pix_str, sprintf(' $%d^2$',2^(l_range(i))*res_ref(1)) ' &'];
    NRMSE_str = [NRMSE_str, sprintf(' %3.2e',NRMSE{end}(i)) ' &'];
    SER_str = [SER_str, sprintf(' %3.1f',10*log10(SER{end}(i))) ' &'];
    maxerror_str = [maxerror_str , sprintf(' %3.2e',maxerror{end}(i)) ' &'];
    time_str = [time_str , sprintf(' %4.1f',t_discrete{end}(i)) ' &'];
    speedup_str = [speedup_str , sprintf(' %4.2f',t_discrete{end}(i)/t_ref(end)) ' &'];
end
fexp3 = fopen(sprintf('../table_exp3_%s.tex',sens_model),'w+');
fprintf(fexp3,'\\begin{center}\n\\begin{tabular}{%s}\n\\hline \\# pixels & %s \\\\\n\\hline\\hline NRMSE & %s \\\\\n\\hline SER (dB) & %s \\\\\n\\hline max.\\ error & %s \\\\\n\\hline time (s) & %s \\\\\n\\hline rel.\\ time & %s \\\\\n\\hline\n\\end{tabular}\n\\end{center}',dcl_str,nb_pix_str(1:end-1),NRMSE_str(1:end-1),SER_str(1:end-1),maxerror_str(1:end-1),time_str(1:end-1),speedup_str(1:end-1));
fclose(fexp3);