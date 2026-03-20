function DATA = lowrankMRF2D_forward(imagesCompressed,Phi,FT,idproj,coilmap,datsize)
% LOWRANKMRF2D_FORWARD Transforms low-rank MRF imagesCompressed back to k-space
%
% Inputs
% 1. imagesCompressed - (double, Ny x Nx x K) the low-rank MRF images,
% where K is the rank of the compressed MRF dictionary.
% 2. Phi -- (double, TRs x K) multiplication with this projection operator transforms data
% from the time domain to the compressed subspace. Multiplication by Phi'
% transforms data from the subspace back to the time domain.
% 3. FT -- NUFFT operator (see NUFFT.m)
% 4. idproj -- (double, TRs x 1) lists the spiral interleaf index for each TR
% 5. coilmap -- (double, Ny x Nx x coils) coil sensitivity maps.
%
% Outputs
% 1. DATA -- (double, read x coils x TRs) MRF k-space data
%
% Jesse Hamilton
% Date: 3/16/2020
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
[N,~,K] = size(imagesCompressed);
nc = size(coilmap,3);
Nex = size(Phi,1);
if nargin < 6
    datsize = [0,0];
end

if isa(FT,'NUFFT')
    imcomp = complex(zeros(N,N,nc,K,'single'));
    for t=1:K
        imcomp(:,:,:,t) = bsxfun(@times,coilmap,imagesCompressed(:,:,t));
    end
    kcomp = FT*imcomp; % k-space data, in compressed subspace
else
    nr = datsize(1);
    np = datsize(2);
    kcomp = complex(zeros(nr,np,nc,K,'single'));
    for t = 1:K
        kcomp(:,:,:,t) = reshape(FT * imagesCompressed(:,:,t),[nr,np,nc]);
    end
end

[nr,np,nc,K]=size(kcomp);
DATA = complex(zeros(nr,nc,Nex,class(kcomp)));
for coil = 1:nc
    ksvd = reshape(kcomp(:,:,coil,:),[nr*np K])*Phi';
    for t=1:Nex
        p = idproj(t);
        DATA(:,coil,t) = ksvd(nr*(p-1)+1:nr*p,t) + squeeze(DATA(:,coil,t));
    end
end
