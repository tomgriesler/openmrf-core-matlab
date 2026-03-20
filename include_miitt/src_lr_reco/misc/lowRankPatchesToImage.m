function [ new_image, s_vals ] = lowRankPatchesToImage(imsize, block_dim, step, shift_idx, U_big,V_big,s_LLR)
%llr_thresh Locally Low Rank regularization
%   Implements Locally Low Rank regularization through singular value soft
%   thresholding. Non-overlapping patches are extracted from the image and
%   reshaped into a matrix. Random shifting is applied before/after the
%   reshaping
%
%  Inputs:
%    alpha [ny, nz, K] -- coefficient images with K coefficients
%    lambda -- soft threshold parameter
%    block_dim [Wy, Wz] -- spatial dimensions of image block
%    randshift -- true to perform random shifting
%
%  Outputs:
%     alpha_thresh [ny, nz, K] -- coefficient image after singular value
%        thresholding
%     s_vals [ny / Wy, nz / Wz, K] -- singular values of each block before
%        thresholding
%
%  Notes:
%     The image dimensions should be divisible by the block sizes
Wy=block_dim(1);
Wx=block_dim(2);
N=imsize(1);
K=imsize(3);
reps=size(s_LLR,2);

alpha_LLR = zeros(Wy*Wx, K, reps);
for ii=1:reps
    UU = U_big(:,:,ii);
    VV = V_big(:,:,ii);
    s2 = s_LLR(:,ii);
    alpha_LLR(:,:,ii) = UU*diag(s2)*VV';
end

counter = zeros(N,N);
i=0;
new_image = zeros(N,N,K);
for lrow=1:step:N-Wy+1
    for lcol=1:step:N-Wx+1
        i=i+1;
        new_image(lrow:lrow+Wy-1,lcol:lcol+Wx-1,:) = new_image(lrow:lrow+Wy-1,lcol:lcol+Wx-1,:) + reshape(alpha_LLR(:,:,i),[Wy Wx K]);
        counter(lrow:lrow+Wy-1,lcol:lcol+Wx-1)=counter(lrow:lrow+Wy-1,lcol:lcol+Wx-1)+1;
    end
end
counter(counter==0)=1;
new_image = bsxfun(@rdivide,new_image,counter);
new_image = circshift(new_image,-shift_idx);


s_vals = permute(reshape(s_LLR, [size(s_LLR,1), sqrt(size(s_LLR,2)), sqrt(size(s_LLR,2))]), [2, 3, 1]);

end
