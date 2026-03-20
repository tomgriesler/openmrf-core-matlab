function [DATA_v, SI_Ratio, weights] = ROVIR(DATA, images_coils, method, numCoils, thresh, orthonormalize, roiSignal, roiInterf)
% Performs ROVir artifact suppression using coil beamforming, as described
% in the following paper:
% Region-optimized virtual (ROVir) coils: Localization and/or suppression of spatial regions using sensor-domain beamforming
% Daeun Kim, et al. MRM 2021.
% https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28706
% 
% code implementation: Jesse Hamilton,  12/20/2021
% code optimization:   Maximilian Gram, 12/20/2025 ... we were both obviously bored before Christmas
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    [NRead, NCoils, NR] = size(DATA);
    Nxy                 = size(images_coils,1);
    
    Xs = reshape(bsxfun(@times,images_coils,roiSignal),[],NCoils);
    Xi = reshape(bsxfun(@times,images_coils,roiInterf),[],NCoils);
    A = Xs' * Xs;
    B = Xi' * Xi;
    clear Xs Xi;
    
    [v,s,~] = eig(A,B);
    svals = abs(diag(s)/s(1));
    if svals(1) < svals(end)
        v = v(:,end:-1:1);
    end
    for i=1:size(v,2)
        v(:,i)=v(:,i)/norm(v(:,i));
    end
    if orthonormalize
        v = gson(v);
    end
    
    retainedSignal=zeros(NCoils,1);
    retainedInterf=zeros(NCoils,1);
    for i=1:NCoils
        retainedSignal(i) = norm(v(:,1:i)*v(:,1:i)'*A*v(:,1:i)*v(:,1:i)','fro')/norm(A,'fro')*100;
        retainedInterf(i) = norm(v(:,1:i)*v(:,1:i)'*B*v(:,1:i)*v(:,1:i)','fro')/norm(B,'fro')*100;
    end
    
    SI_Ratio = retainedSignal./retainedInterf;
    
    if strcmp(method, 'manual')
        keeps     = 1:numCoils;
        eigenvecs = v(:,keeps);
    
    elseif strcmp(method, 'auto')
        weights  = repmat(v,[1,1,Nxy,Nxy]);
        weights  = permute(weights,[3,4,1,2]);
        imcoil_v = squeeze(sum(bsxfun(@times,images_coils,weights),3));
        for u=1:NCoils
            interf(u) = sum(sum(abs(bsxfun(@times,imcoil_v(:,:,u),roiInterf).^2)));
            sig(u)    = sum(sum(abs(bsxfun(@times,imcoil_v(:,:,u),roiSignal).^2)));
        end
        sigratio  = sig./interf;
        keeps     = find(sigratio >= thresh);
        eigenvecs = v(:,keeps);
    end
    
    DATA_v = zeros(NRead, size(eigenvecs,2), NR, 'like', DATA);
    for r = 1:NR
        DATA_v(:,:,r) = DATA(:,:,r) * eigenvecs;
    end
    
    disp(['ROVir: keep ' num2str(size(DATA_v,2)) ' of ' num2str(size(DATA,2)) ' coils' ]);

end