clear all;clc;close all;
% Loading precomputed experiment data
load brain_recon_data.mat;

% preparing sensitivity data for reconstruction
res = size(reference);
NbCoils = length(sens);
[X,Y] = GenerateFullCart2DGrid(res);
R = [X(support)'/res(1);Y(support)'/res(2)];clear X Y;
M = Sinusoidal2DMatrix(R,sens{1}.param);clear R;
param.sensitivity = zeros([res,NbCoils]);
SoS = zeros(res);
for c=1:NbCoils
    tmp = zeros(res);
    tmp(support) = M*sens{c}.data;
    param.sensitivity(:,:,c) = tmp;
    SoS = SoS + abs(tmp).^2;
end
SoS(~support)=1;
SoS = ones(res);
clear M tmp c;
x_init = zeros(res);%reference;

% Setting reconstruction parameters
tol = 1e-7;
NitStep = 50;

lambda_opt = zeros(2,3,3,2);
% CG
lambda_opt(1,:,1,1) = [7 6 6]*1e-4; % spiral
lambda_opt(1,:,2,1) = [1.1 1.2 1.2]*1e-4;
lambda_opt(1,:,3,1) = [3 3.5 4]*1e-5;
lambda_opt(2,:,1,1) = [5 7 8]*1e-5; % EPI
lambda_opt(2,:,2,1) = [5 7 6]*1e-6;
lambda_opt(2,:,3,1) = [0 5 0]*1e-7;
% TV
lambda_opt(1,:,1,2) = [4 4 3]*1e-4; % spiral
lambda_opt(1,:,2,2) = [8 7 8]*1e-5;
lambda_opt(1,:,3,2) = [3 2 3]*1e-5;
lambda_opt(2,:,1,2) = [1.5 1 1.5]*1e-4; % EPI
lambda_opt(2,:,2,2) = [3.5 3 4]*1e-5;
lambda_opt(2,:,3,2) = [1.5 1 2]*1e-5;

FD = @(x) {x-circshift(x,[1,0]),x-circshift(x,[0,1])};
FDT = @(g) g{1}-circshift(g{1},[-1,0])+g{2}-circshift(g{2},[0,-1]);

% performing reconstructions
SER_opt = zeros(2,3,3,2);
for id_traj = 1:length(traj_list)
    traj = traj_list{id_traj};
    switch traj
        case 'spiral'
            D = ones(res);
        case 'EPI'
            D = 1./sqrt(SoS);
    end
    Pre = @(vx) reshape(vx,res).*D;
    Post = @(x) reshape(x.*D,prod(res),1);
    param.G = PSF{id_traj};
    ATA = @(x) support.*EHE(support.*x,param);
    for id_simu = 1:length(simu_list)
        simu = simu_list{id_simu};
        for id_SNR = 1:length(SNR_list)
            SNR = SNR_list(id_SNR);
            %% Tikhnov (CG)
            M = @(x) ATA(x)+lambda_opt(id_traj,id_simu,id_SNR,1)*x;
            PMP = @(vx) Post(M(Pre(vx)));
            flag = 1;
            counter = 0;
            fprintf('%s %s %d dB. %s lbd=%.2e\n',traj,simu,SNR,'CG',lambda_opt(id_traj,id_simu,id_SNR,1));
            vx_init = Post(x_init./(D).^2);
            va = Post(adjoint{id_traj,id_simu,id_SNR});
            nva = norm(va);
            resvectot = [];
            while flag
                [vx,flag,relres,iter,resvec] = pcg(PMP,va,tol,NitStep,[],[],vx_init);
                x = Pre(vx);
                if iter==0
                    break;
                end
                resvectot = [resvectot; resvec(2:end)/nva];
                counter = counter+iter;
                SER = 20*log10(norm(reference(:),2)/norm(reference(:)-x(:),2));
                fprintf('iteration %3d. relres %5.3e (%d). SER %5.3f\n',counter,relres,flag,SER);
                semilogy(resvectot);drawnow;
                vx_init = vx;
            end
            SER_opt(id_traj,id_simu,id_SNR,1) = SER;
            save(sprintf('Brain_%s_%s_%ddB_%s_%d-%ddB.mat',traj,simu,SNR,'CG',floor(SER),round(100*(SER-floor(SER)))),'x','SER');
            x_print = uint16((2^16-1)*x);
            imwrite(x_print,sprintf('Brain_%s_%s_%ddB_%s_%d-%ddB.png',traj,simu,SNR,'CG',floor(SER),round(100*(SER-floor(SER)))));
            %% TV (IRLS)
            epsilon = 3e-7;
            counter = 0;
            fprintf('%s %s %d dB. %s lbd=%.2e\n',traj,simu,SNR,'TV',lambda_opt(id_traj,id_simu,id_SNR,2));
            for it = 1:15
                FDx = FD(x);
                w = abs(FDx{1}).^2+abs(FDx{2}).^2;
                w = 1./sqrt(w+epsilon*max(w(:)));
                WFD = @(w,x) {w.*(x-circshift(x,[1,0])),w.*(x-circshift(x,[0,1]))};
                M = @(x) ATA(x)+lambda_opt(id_traj,id_simu,id_SNR,2)*FDT(WFD(w,x))/2;
                PMP = @(vx) Post(M(Pre(vx)));
                [vx,flag,relres,iter,resvec] = pcg(PMP,va,tol,30,[],[],vx);
                x = Pre(vx);
                if iter==0
                    break;
                end
                counter = counter+iter;
                SER = 20*log10(norm(reference(:),2)/norm(reference(:)-x(:),2));
                fprintf('iteration %3d. relres %5.3e (%d). SER %5.3f\n',counter,relres,flag,SER);
                semilogy(resvectot);drawnow;
            end
            SER_opt(id_traj,id_simu,id_SNR,2) = SER;
            save(sprintf('Brain_%s_%s_%ddB_%s_%d-%ddB.mat',traj,simu,SNR,'TV',floor(SER),round(100*(SER-floor(SER)))),'x','SER');
            x_print = uint16((2^16-1)*x);
            imwrite(x_print,sprintf('Brain_%s_%s_%ddB_%s_%d-%ddB.png',traj,simu,SNR,'TV',floor(SER),round(100*(SER-floor(SER)))));
        end
    end
end

%% Writting error maps
for id_traj = 1:length(traj_list)
    traj = traj_list{id_traj};
    for id_simu = 1:length(simu_list)
        simu = simu_list{id_simu};
        for id_SNR = 1:length(SNR_list)
            SNR = SNR_list(id_SNR);
            %% Tikhnov (CG)
            load(sprintf('Brain_%s_%s_%ddB_%s_%d-%ddB.mat',traj,simu,SNR,'CG',floor(SER_opt(id_traj,id_simu,id_SNR,1)),round(100*(SER_opt(id_traj,id_simu,id_SNR,1)-floor(SER_opt(id_traj,id_simu,id_SNR,1))))));
            x_print = uint16((2^16-1)*(abs(x-reference)/m).^gamma);
            m = [m, max(abs(x(:)-reference(:)))];
            %% TV (IRLS)
            SER = SER_opt(id_traj,id_simu,id_SNR,2);
            load(sprintf('Brain_%s_%s_%ddB_%s_%d-%ddB.mat',traj,simu,SNR,'TV',floor(SER),round(100*(SER-floor(SER)))));
            m = [m, max(abs(x(:)-reference(:)))];
        end
    end
end
m = max(m);
gamma = 0.8;
%m = 0.5523;
for id_traj = 1:length(traj_list)
    traj = traj_list{id_traj};
    for id_simu = 1:length(simu_list)
        simu = simu_list{id_simu};
        for id_SNR = 1:length(SNR_list)
            SNR = SNR_list(id_SNR);
            %% Tikhnov (CG)
            load(sprintf('Brain_%s_%s_%ddB_%s_%d-%ddB.mat',traj,simu,SNR,'CG',floor(SER_opt(id_traj,id_simu,id_SNR,1)),round(100*(SER_opt(id_traj,id_simu,id_SNR,1)-floor(SER_opt(id_traj,id_simu,id_SNR,1))))));
            x_print = uint16((2^16-1)*(abs(x-reference)/m).^gamma);
            %m = [m, max(abs(x(:)-reference(:)))];
            imwrite(x_print,sprintf('Brain_%s_%s_%ddB_%s_%d-%ddB_Error.png',traj,simu,SNR,'CG',floor(SER_opt(id_traj,id_simu,id_SNR,1)),round(100*(SER_opt(id_traj,id_simu,id_SNR,1)-floor(SER_opt(id_traj,id_simu,id_SNR,1))))));
            %% TV (IRLS)
            SER = SER_opt(id_traj,id_simu,id_SNR,2);
            load(sprintf('Brain_%s_%s_%ddB_%s_%d-%ddB.mat',traj,simu,SNR,'TV',floor(SER),round(100*(SER-floor(SER)))));
            x_print = uint16((2^16-1)*(abs(x-reference)/m).^gamma);
            %m = [m, max(abs(x(:)-reference(:)))];
            imwrite(x_print,sprintf('Brain_%s_%s_%ddB_%s_%d-%ddB_Error.png',traj,simu,SNR,'TV',floor(SER),round(100*(SER-floor(SER)))));
        end
    end
end