function imcomp = lowrankMRF2D_adjoint(DATA,coilmap,idproj,Phi,FT,np)
% LOWRANKMRF2D_ADJOINT Compresses MRF k-space along the time dimension,
% then grids the compressed k-space to the image domain. Yields "low-rank"
% images.
% 
% Input
% 1. DATA -- (double, read x coils x TRs) the MRF k-space data
% 2. coilmap -- (double, Ny x Nx x coils) coil sensitivity maps
% 3. idproj -- (double, TRs x 1) spiral interleaf index for each TR
% 4. Phi -- (double, TRs x K) projection matrix that converts data from the
%               time domain to the compressed subspace. K is the rank of
%               the dictionary after compression.
% 5. FT -- NUFFT operator (see NUFFT.m)
% 6. np -- (double, 1x1) the number of spiral interleaves in the fully
% sampled data
% 
% Outputs
% imcomp -- (double, Ny x Nx x K) the low-rank images (i.e., the images in
% the compressed subspace).
% 
% Jesse Hamilton
% Date: 3/16/2020
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

[nr,nc,Nex] = size(DATA);
N = size(coilmap,1);
K = size(Phi,2);

% Compress k-space along the time dimension
ksvd = complex(zeros(nr*np,Nex,class(DATA)));
kcomp = complex(zeros(nr,np,nc,K,class(DATA)));
for coil = 1:nc
    for t=1:Nex
        p = idproj(t); % spiral interleaf
        ksvd(nr*(p-1)+1:nr*p,t) = DATA(:,coil,t);
    end
    kcomp(:,:,coil,:) = reshape(ksvd*Phi,[nr np K]);
end

% Grid the compressed k-space data to the image domain to yield "low-rank
% images"
imcomp = complex(zeros(N,N,K,class(DATA)));
if isa(FT,'NUFFT')
    for t=1:K
        imcomp(:,:,t) = sum(conj(coilmap).*(FT'*kcomp(:,:,:,t)),3);%./sum(abs(coilmap).^2,3);
    end
else
    for t=1:K
        imcomp(:,:,t) = FT' * reshape(kcomp(:,:,:,t),[nr*np nc]);
    end
end
