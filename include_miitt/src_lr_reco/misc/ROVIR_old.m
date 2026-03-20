function [DATA_v, SI_Ratio, weights] = ROVIR_old(DATA, imcoil, N, os,  thresh, method, ROIType, orthonormalize,roiSignal,roiInterf)
% Performs ROVir artifact suppression using coil beamforming, as described
% in the following paper:
% Region-optimized virtual (ROVir) coils: Localization and/or suppression of spatial regions using sensor-domain beamforming
% Daeun Kim, et al. MRM 2021.
% https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28706
% 
% code implementation: Jesse Hamilton, 12/20/2021
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% DATA -- kspace data, size [readout x coils x TRs]
% idproj -- projection index for each TR

if ndims(DATA)==3
    [nr,nc,Nex] = size(DATA);
else
    [nr,np,nc,Nex] = size(DATA);
end
SI_Ratio=[];


%====== automatic outer FOV suppression ========== %
% % %         % circular support
% buffer = 0.1;
% [X,Y] = meshgrid( linspace(-1,1,N*os), linspace(-1,1,N*os) );
% distance = sqrt(X.^2 + Y.^2);
% roiSignal = zeros(N*os,N*os);
% roiSignal(distance < 0.5)=1;
% roiInterf = zeros(N*os,N*os);
% roiInterf(distance > 0.5+buffer) = 1;

imSOS = makesos(imcoil,3);
% %     % %     % square support
if strcmp(ROIType,'auto')
%     buffer=1;
%       buffer=8;
    buffer = ceil(N*0.05);
    
% % %     % CIRCULAR REGION OF SUPPORT
%     [X,Y] = meshgrid( linspace(-1,1,N*os), linspace(-1,1,N*os) );
%     distance = sqrt(X.^2 + Y.^2);
%     roiSignal = zeros(N*os,N*os);
%     roiSignal(distance < 0.5*sqrt(2))=1;
% % %     roiSignal(distance < 0.5*sqrt(2))=1;
%     roiInterf = zeros(N*os,N*os);
%     roiInterf(distance>0.5*sqrt(2)+.05) = 1;
% %     roiInterf(distance > 0.5+buffer/N) = 1;
% 
% %     % SQUARE REGION OF SUPPORT
    roiSignal = zeros(N*os,N*os);
    roiSignal(N*os/2-N/2+1:N*os/2+N/2,N*os/2-N/2+1:N*os/2+N/2)=1;
    roiInterf = zeros(N*os,N*os);
    roiInterf(N*os/2-N/2+1-buffer:N*os/2+N/2+buffer,N*os/2-N/2+1-buffer:N*os/2+N/2+buffer) = 1;
    roiInterf = 1-roiInterf;
%     
    
    
    
    
    
%     
%     imInterf = imSOS.*roiInterf;
%     meanValue = mean(abs(imSOS(roiSignal==1)))
%     idx = find(abs(imSOS.*roiInterf) > .4*meanValue);
%     imageInterf = 0*imInterf;
%     imageInterf(idx) = 1;
% %     figure,imshow(imageInterf,[])
%     roiInterf = imageInterf;
    
%     figure,subplot(131);imshow(imSOS.^.5.*roiSignal,[]);
%     subplot(132); imshow(imSOS.^.5.*roiInterf,[]);
%     subplot(133); imshow(imSOS.^.5.*(roiSignal+roiInterf).^.5,[]);
%     drawnow
    
    % %====== manual outer FOV suppression ========%
% else
%     figure,imagesc(makesos(imcoil,3).^.5);axis image; axis off; colormap gray; title('draw ROI over signal area');
%     h = imellipse;
%     pause;
%     roiSignal = h.createMask;
% %     roiSignal = roipoly;
% 
%     if 1
% %         figure,imagesc(makesos(imcoil,3).^.5);axis image; axis off; colormap gray; 
%         title('draw ROI over interference area #1');
%         h2 = imellipse;
%         pause;
%         roiInterf = 1-h2.createMask;
%     
% % %         figure,imagesc(makesos(imcoil,3).^.5);axis image; axis off; colormap gray; 
% %         title('draw ROI over interference area #2');
% %         h3 = imellipse;
% %         pause;
% %         roiInterf2 = h3.createMask;
% % 
% %         roiInterf = 0*roiInterf1;
% %         roiInterf(roiInterf1)=1;
% %         roiInterf(roiInterf2)=1;
%     
%     else
%         roiInterf = 1-imdilate(roiSignal,strel('disk',round(N*0.1)));
%         % % roiInterf = roipoly;
%         %     roiInterf = 1-imdilate(roiSignal,strel('disk',10));
%     
%         imInterf = imSOS.*roiInterf;
%         meanValue = mean(abs(imSOS(roiSignal==1)))
%         idx = find(abs(imSOS.*roiInterf) > .4*meanValue);
%         imageInterf = 0*imInterf;
%         imageInterf(idx) = 1;
%         %     figure,imshow(imageInterf,[])
%         roiInterf = imageInterf;
%     end
    
    % %         title('click on center of Signal Region');
    % %         radius = 40;
    % %         buffer = 10;
    % %         [xi,yi]=ginput(1);
    % %         xi=round(xi); yi=round(yi);
    % %         [X,Y]=meshgrid(1:N*os,1:N*os);
    % %         distance = sqrt((X-xi).^2 + (Y-yi).^2);
    % %         roiSignal = zeros(N*os,N*os);
    % %         roiSignal(distance < radius)=1;
    % %         roiInterf = zeros(N*os,N*os);
    % %         roiInterf(distance > radius+buffer) = 1;
%     figure,subplot(131);imshow(imSOS.^.5.*roiSignal,[]);
%     subplot(132); imshow(imSOS.^.5.*roiInterf,[]);
%     subplot(133); imshow(imSOS.^.5.*(roiSignal+roiInterf).^.5,[]);
%     drawnow
end

% g = reshape(permute(imcoil,[3 1 2]),nc,[]);

% % % % % ORIGINAL WAY
% idx_signal = find(roiSignal==1);
% A = conj(g(:,idx_signal))*g(:,idx_signal).';
% idx_interf = find(roiInterf==1);
% B = conj(g(:,idx_interf))*g(:,idx_interf).';

Xs = reshape(bsxfun(@times,imcoil,roiSignal),[],nc);
Xi = reshape(bsxfun(@times,imcoil,roiInterf),[],nc);
A = Xs' * Xs;
B = Xi' * Xi;

[v,s,d] = eig(A,B);
svals = abs(diag(s)/s(1));
if svals(1) < svals(end)
    v = v(:,end:-1:1);
    svals = svals(end:-1:1);
end
%     svals = svals/max(svals)

% vnorm = v;
% vnorm = gson(v);

for i=1:size(v,2)
    v(:,i)=v(:,i)/norm(v(:,i));
end
if orthonormalize
    v = gson(v);
end
% vnorm = gsog(v);


retainedSignal=zeros(nc,1);
retainedInterf=zeros(nc,1);
for i=1:nc
    retainedSignal(i) = norm(v(:,1:i)*v(:,1:i)'*A*v(:,1:i)*v(:,1:i)','fro')/norm(A,'fro')*100;
    retainedInterf(i) = norm(v(:,1:i)*v(:,1:i)'*B*v(:,1:i)*v(:,1:i)','fro')/norm(B,'fro')*100;
end

% retainedSignal = retainedSignal/max(retainedSignal);
% retainedInterf = retainedInterf/max(retainedInterf);
SI_Ratio = retainedSignal./retainedInterf;

% figure; plot(retainedSignal,'linewidth',3); hold all; plot(retainedInterf,'linewidth',3);
% plot(SIR,'linewidth',3);
% legend('Signal','Interf','SIR');
% drawnow


% fprintf('\n');
% fprintf('Signal Retained: %.1f%%\n',retainedSignal(nv)*100);
% fprintf('Interference Retained: %.1f%%\n',retainedInterf(nv)*100);


% signalEnergy = zeros(nc,1);
% interfEnergy = zeros(nc,1);
% parfor i=1:nc
%     weights = v(:,1:i);
%     weights = repmat(weights,[1,1,N*os,N*os]);
%     weights = permute(weights,[3,4,1,2]);
%     imVirtualCoils = squeeze(sum(bsxfun(@times,imcoil,weights),3)); % Y, X, Coils
% %     imVirtualCoils = permute(imVirtualCoils,[3 1 2]); % coils, Y, X
% %     imcoilSignal = imVirtualCoils(:,idx_signal); % ORIGINAL WAY
% %     imcoilInterf = imVirtualCoils(:,idx_interf);
%     imcoilSignal = reshape(bsxfun(@times,imVirtualCoils,roiSignal),N*N*os*os,[]);
%     imcoilInterf = reshape(bsxfun(@times,imVirtualCoils,roiInterf),N*N*os*os,[]);
% 
%     signalEnergy(i) = sum(sum(abs(imcoilSignal)));
%     interfEnergy(i) = sum(sum(abs(imcoilInterf)));
% %     signalEnergy(i) = mean(mean(abs(imcoilSignal)));
% %     interfEnergy(i) = mean(mean(abs(imcoilInterf)));
% %     if i < nc
% %         weights = v(:,i+1:end);
% %         weights = repmat(weights,[1,1,N*os,N*os]);
% %         weights = permute(weights,[3,4,1,2]);
% %         imcoilInterf = squeeze(sum(bsxfun(@times,imcoil,weights),3));
% %         imcoilInterf = permute(imcoilInterf,[3 1 2]);
% %         imcoilInterf = imcoilInterf(:,idx_interf);
% %         interfEnergy(i) = sum(sum(abs(imcoilInterf).^2));
% %     else
% %         weights = v;
% %         weights = repmat(weights,[1,1,N*os,N*os]);
% %         weights = permute(weights,[3,4,1,2]);
% %         imcoilInterf = squeeze(sum(bsxfun(@times,imcoil,weights),3));
% %         imcoilInterf = permute(imcoilInterf,[3 1 2]);
% %         imcoilInterf = imcoilInterf(:,idx_interf);
% %         interfEnergy(i) = sum(sum(abs(imcoilInterf).^2));
% %     end
% end
% signalEnergy = signalEnergy/max(signalEnergy)*100;
% interfEnergy = interfEnergy/max(interfEnergy)*100;
% SI_Ratio = signalEnergy./interfEnergy;
% % % SI_Ratio = SI_Ratio/min(SI_Ratio);
% 

% disp('Coils   Signal   Interf   S/I');
% disp([(1:nc)',signalEnergy(:), interfEnergy(:),SI_Ratio(:)])

if strcmp(method,'numCoils')
    nv = thresh;
    keeps = 1:nv;
    dels = nv+1:nc;
    eigenvecs = v(:,keeps);
    eigenvecs_i = v(:,dels);
% elseif strcmp(method,'interference')
%     nv = find(retainedInterf <= thresh, 1, 'last');
% elseif strcmp(method,'signal')
%     nv = find(retainedSignal <= thresh, 1, 'last');
% elseif strcmp(method,'SIRatio')
%     nv = find(SI_Ratio >= thresh, 1, 'last');
%     while nv <= 1
%         thresh = thresh / 2;
%         nv = find(SI_Ratio >= thresh, 1, 'last');
%     end
% elseif strcmp(method,'svals')
%     svals = svals/svals(1);
% %     svalsEnergy = cumsum(svals.^2./sum(svals.^2));
% %     K = find(svalsEnergy > 0.9999,1,'first');
%     nv = find(svals<thresh,1,'first');
%     if isempty(nv)
%         nv=nc;
%     end
% elseif strcmp(method,'svalsEnergy')
%     svals = svals/svals(1);
%     svalsEnergy = cumsum(svals.^2./sum(svals.^2));
%     nv = find(svalsEnergy > thresh,1,'first');
%     if isempty(nv)
%         nv=nc;
%     end
elseif strcmp(method,'auto-thresh')
    weights = repmat(v,[1,1,N*os,N*os]);
    weights = permute(weights,[3,4,1,2]);
    imcoil_v = squeeze(sum(bsxfun(@times,imcoil,weights),3));
    for u=1:nc
        interf(u) = sum(sum(abs(bsxfun(@times,imcoil_v(:,:,u),roiInterf).^2)));
        sig(u)    = sum(sum(abs(bsxfun(@times,imcoil_v(:,:,u),roiSignal).^2)));
    end
    sigratio = sig./interf;
    disp([sig(:) interf(:) sigratio(:)])
    % disp(interf(:)/max(interf(:)))
    keeps = find(sigratio >= thresh);
    dels = find(sigratio < thresh);
    nv = length(keeps);
    eigenvecs = v(:,keeps);
    eigenvecs_i = v(:,dels);

    imcoil_v=mat2gray(abs(imcoil_v));

    figure,
    subplot(121); display3D(abs(bsxfun(@times,imcoil_v,roiSignal)).^.5,ceil(sqrt(nc)),ceil(sqrt(nc)),[0 1],'doNorm',0); 
    subplot(122); display3D(abs(bsxfun(@times,imcoil_v,roiInterf)).^.5,ceil(sqrt(nc)),ceil(sqrt(nc)),[0 1],'doNorm',0); 
    colormap(gray)
    drawnow;

    fprintf('Removing coils: \n')
    disp(dels)
end
fprintf('KEEPING %d VIRTUAL COILS\n',nv);



weights = repmat(eigenvecs,[1,1,N*os,N*os]);
%         weights = repmat(v,[1,1,N*os,N*os]);
weights = permute(weights,[3,4,1,2]);
imcoil_v = squeeze(sum(bsxfun(@times,imcoil,weights),3));


%         weights = repmat(eigenvecs,[1,1,N*os,N*os]);
%         weights = repmat(v,[1,1,N*os,N*os]);
%         weights = permute(weights,[3,4,1,2]);
weights_i = repmat(eigenvecs_i,[1,1,N*os,N*os]);
%     weights_i = repmat(v(:,end-nv+1:end),[1,1,N*os,N*os]);
weights_i = permute(weights_i,[3,4,1,2]);
imcoil_i = squeeze(sum(bsxfun(@times,imcoil,weights_i),3));



if ndims(DATA)==3
    weights = repmat(eigenvecs,[1,1,nr,Nex]);
    weights = permute(weights,[3,1,4,2]);
    DATA_v = squeeze(sum(bsxfun(@times,DATA,weights),2));
    DATA_v = permute(DATA_v,[1,3,2]);
elseif ndims(DATA)==4
    weights = repmat(eigenvecs,[1,1,nr,np,Nex]);
    weights = permute(weights,[3,4,1,5,2]);
    DATA_v = squeeze(sum(bsxfun(@times,DATA,weights),3)); % read proj meas coils
    DATA_v = reshape(DATA_v,[nr,np,Nex,nv]);
    DATA_v = permute(DATA_v,[1,2,4,3]); % read proj coils meas
else
    error('kspace data has wrong dimensions')
end

maxval = max(max(makesos(imcoil,3).^.75));
figure
subplot(234); imshow(makesos(imcoil,3).^.75,[]); title('All Coils');
subplot(232); imshow(roiSignal.*makesos(imcoil,3).^.75,[]); title('All Coils (signal mask)');
subplot(233); imshow(roiInterf.*makesos(imcoil,3).^.75,[]); title('All Coils (interf mask)');
subplot(235); imshow(makesos(imcoil_v,3).^.75,[]); title('ROVir - Signal');
subplot(236); imshow(makesos(imcoil_i,3).^.75,[]); title('ROVir - Interf');
drawnow