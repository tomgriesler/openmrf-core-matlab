function [s_LLR,U_big,V_big,shift_idx] = imageToLowRankPatches(alpha, block_dim,step, shift_idx)

[N,~,K]=size(alpha);
Wy = block_dim(1);
Wx = block_dim(2);
% Wy = 6;
% Wx = 6; % block size

reps=length(1:step:N-Wy+1)^2;

if nargin < 4
shift_idx = [randi(Wy) randi(Wx) 0];
end
alpha = circshift(alpha,shift_idx);

alpha_LLR = complex(zeros(Wy*Wx,K,reps,class(alpha)));
i=0;
for lrow=1:step:N-Wy+1
    for lcol=1:step:N-Wx+1
        i=i+1;
        alpha_LLR(:,:,i)=reshape(alpha(lrow:lrow+Wy-1,lcol:lcol+Wx-1,:),[Wy*Wx K]);
    end
end


% alpha_LLR = permute(alpha_LLR, [1 3 2]);
s_LLR = zeros(K, reps);
U_big = zeros(Wy*Wx, K, reps);
V_big = zeros(K,K,reps);

if Wy*Wx < K
    s_LLR = zeros(Wy*Wx,reps);
    U_big = zeros(Wy*Wx,Wy*Wx,reps);
    V_big = zeros(K,Wy*Wx,reps);
end

% threshold singular values

for ii=1:reps
    [UU, SS, VV] = svd(alpha_LLR(:,:,ii), 'econ');
    s_LLR(:,ii) = diag(SS);
    %     s2 = SoftThresh(s_LLR(:,ii), lambda);
    %     alpha_LLR(:,:,ii) = UU*diag(s2)*VV';
    U_big(:,:,ii) = UU;
    V_big(:,:,ii) = VV;
end