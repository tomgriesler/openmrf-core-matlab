function [x,t0,dx,obj,update] = nonlinearCGDescent(x0,dx,E,Et,y,params)

time0=clock;
[ny,nx,K] = size(x0);

maxlsiter = 6;
gradToll = 1e-8;
if ~isfield(params,'betaMethod')
    betaMethod = 'Dai-Yuan';
else
    betaMethod = params.betaMethod;
end
if ~isfield(params,'maxDisplay')
    maxDisplay = 1;
else
    maxDisplay = params.maxDisplay;
end
if ~isfield(params,'beta')
    params.beta = 0.8;
end
if ~isfield(params,'t0')
    params.t0=1;
end
if ~isfield(params,'frameidx')
    params.frameidx = 1:K;
end
if ~isfield(params,'stopthresh')
    params.stopthresh=Inf;
end
if ~isfield(params,'alpha')
    params.alpha = 0.01;
end
if ~isfield(params,'maxShowValue')
    params.maxShowValue = 0.85;
end
frameidx=params.frameidx;
beta = params.beta;
t0 = params.t0;
alpha= params.alpha;

x = x0;
update = Inf;

shift_idx = [randi(params.block_dim(1)) randi(params.block_dim(1)) 0];

g0 = grad(x,y,E,Et,params,shift_idx);
if isempty(dx)
    dx = -g0;
end

iter = 0;
obj = [];

% figure

while(1)
    
    f0 = objective(x,dx,0,y,E,params,shift_idx);
    t = t0;
    f1 = objective(x,dx,t,y,E,params,shift_idx);
    lsiter = 0;
    while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)
        lsiter = lsiter + 1;
        t = t * beta;
        f1 = objective(x,dx,t,y,E,params,shift_idx);
    end
    
    % control the number of line searches by adapting the initial step search
    if lsiter > 2, t0 = t0 * beta;end
    if lsiter < 1, t0 = t0 / beta; end
    
    xprev=x;
    x = (x + t*dx);
    
    % print some numbers for debug purposes
    % update=norm1(x(:) - xprev(:))/numel(x(:));

    if iter > 0
        update = (obj(iter) - f1) / obj(iter);
    end

    %     if mod(iter+1,params.verbose)==0
    fprintf('%d,    obj %g, L-S %d, t0 %f, norm %g\n',iter,f1,lsiter,t0,update);
    displayElapsedTime(time0,clock);
    %     end
    
    %     if t0 < stopThresh, break; end
    if abs(update) < params.stopThresh
        disp('met stopping criterion - quitting')
        break;
    end
    
    if mod(iter+1,params.updateFreq)==0
        numframes = length(frameidx);
        % %         frameidx = round(linspace(1,K,numframes));
        %         imbig=zeros(length(params.pixelRangeToShow),length(params.pixelRangeToShow),numframes);
        %         for i=1:numframes
        %             imbig(:,:,i) = mat2gray(abs(x(params.pixelRangeToShow,params.pixelRangeToShow,frameidx(i))));
        %         end
        %         figure(100);clf;display_3D_image(imbig,ceil(sqrt(numframes)),ceil(sqrt(numframes)),[0 maxDisplay],'doNorm',true); colormap(gray(256));
        %         title(sprintf('ite %d: %.4f',iter+1,update)); drawnow
        
%         if isempty(params.Phi)
        figure(100);clf;
        display3D(abs(x(params.pixelRangeToShow,params.pixelRangeToShow,params.frameidx)),ceil(sqrt(numframes)),ceil(sqrt(numframes)),[0 params.maxShowValue]); colormap(gray);
        title(sprintf('Iteration %d / %d: delta=%.4f',iter,params.numIter,update));
        drawnow;
%         else
%             imrec = reshape((uk*reshape(x0,[N*os*N*os,rankk]).').',[N*os N*os nframes]);
%             figure(100);clf;
%             display3D(abs(imrec(params.pixelRangeToShow,params.pixelRangeToShow,params.frameidx)).^.8,ceil(sqrt(numframes)),ceil(sqrt(numframes)),[0 1]); colormap(gray);
%             title(sprintf('Iteration %d / %d: delta=%.4f',iter,params.numIter,update));
%             drawnow;
%             clear imrec
%         end
    end
    
    iter = iter + 1;
    %     obj(iter) = feval(calc_f,x(:));
    %     obj(iter) = calc_F_2DMRI(x(:),y(:),E,params.lambdaSpatialTV,params.lambdaWav,size(x));
    %     obj(iter) = calc_L2part_2DMRF(x(:),y(:),E,size(x));
    obj(iter) = f1;
    
    if (iter >= params.numIter) || (norm(dx(:)) < gradToll), break;end
    if abs(update) <= params.stopThresh
        disp('stopping criteria met')
        break
    end
    
    %conjugate gradient calculation
    g1 = grad(x,y,E,Et,params,shift_idx);
    
    switch betaMethod
        case 'Hestenes-Stiefel'
            bk = g1(:)'*(g1(:)-g0(:)) / (dx(:)'*(g1(:)-g0(:))); % Hestenes-Stiefel
        case 'Dai-Yuan'
            bk = g1(:)'*g1(:) / (dx(:)'*(g1(:)-g0(:))); % Dai-Yuan
        case 'Fletcher-Reeves'
            bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps); % Fletcher-Reeves (default)
        case 'Polak-Ribiere'
            bk = g1(:)'*(g1(:)-g0(:))/(g0(:)'*g0(:)+eps); % Polak-Ribiere
            %             bk= max(bk,0);
        case 'SteepestDescent'
            bk=0;
    end
    
    restartScore = abs(g1(:)'*g0(:)) / abs(g0(:)'*g0(:));
    if restartScore < 0.1
        %         disp('steepest descent')
        bk=0;
    end
    
    g0 = g1;
    dx =  - g1 + bk* dx;
    shift_idx = [randi(params.block_dim(1)) randi(params.block_dim(1)) 0];
end
end

function res = objective(x,dx,t,y,E,params,shift_idx) %**********************************

l1Smooth = 1e-15;
[ny,nx,K] = size(x);
TVOPY = TV_Temp_AnyDim(1);
TVOPX = TV_Temp_AnyDim(2);
TVOPT = TV_Temp_AnyDim(3);
% FTOPT = FT_Temp_AnyDim(3);
ny_diad = 2^nextpow2(ny);

% L2-norm part
w=E((x+t*dx))-y;
% L2Obj = sum(reshape(conj(w).*w,[],1));
L2Obj = sum(vec(conj(w).*w));
% L2Obj = sum(sum(sum(sum(conj(w).*w))));
% L2Obj=w(:)'*w(:);

% total variation
TVObjX = 0;
if params.lambdaSpatialTV > 0
    w = TVOPX * (x + t*dx);
    TVObjX = sum((w(:).*conj(w(:)) + l1Smooth).^(1/2));
end

TVObjY = 0;
if params.lambdaSpatialTV > 0
    w = TVOPY * (x + t*dx);
    TVObjY = sum((w(:).*conj(w(:)) + l1Smooth).^(1/2));
end

TVObjT = 0;
if params.lambdaTemporalTV > 0
    if ~isempty(params.Phi)
        Nex=size(params.Phi,1);
        I = reshape(reshape(x+t*dx,[ny*nx,K])*params.Phi',[ny,nx,Nex]);
        w = TVOPT * I;
        w = reshape(reshape(w,[ny*nx,Nex])*params.Phi,[ny,nx,K]);
        %     w = TVOPT * (x + t*dx);
        TVObjT = sum((w(:).*conj(w(:)) + l1Smooth).^(1/2));
    else
        w = TVOPT*(x+t*dx);
        TVObjT = sum((w(:).*conj(w(:)) + l1Smooth).^(1/2));
    end
end

FTObjT = 0;
if params.lambdaTemporalFT > 0
    w = FTOPT * (x + t*dx);
    FTObjT = sum((w(:).*conj(w(:)) + l1Smooth).^(1/2));
end

% Wavelet
WavObj_dim1 = zeros(K,1);
if params.lambdaWav
    parfor meas=1:K
        WavOP=Waveletxyt('Daubechies',4,6,[ny,nx]);
        w = WavOP*zpad(x(:,:,meas)+t*dx(:,:,meas),[ny_diad ny_diad]);
        %         w=WavOP*(x(:,:,meas)+t*dx(:,:,meas));
        WavObj_dim1(meas) = sum((w(:).*conj(w(:))+l1Smooth).^(1/2));
    end
end
WavObj_dim1 = sum(WavObj_dim1);

% Local low rank image patches
llrObj = 0;
if params.lambdaLLR
    NN = ceil(ny / params.block_dim(1)) * params.block_dim(1);
    w = imageToLowRankPatches(zpad(x+t*dx,NN,NN,K), params.block_dim, params.block_step,shift_idx);
    %     w = my_image_to_llr(x+t*dx,params.block_dim,params.block_step,shift_idx);
    %     [w] = image_to_llr(x+t*dx,params.block_dim,1);
    llrObj = sum((w(:).*conj(w(:))+l1Smooth).^(1/2));
end

res = L2Obj + params.lambdaWav*WavObj_dim1 + ...
    params.lambdaSpatialTV*TVObjX +  params.lambdaSpatialTV*TVObjY + ...
    params.lambdaTemporalTV*TVObjT +  params.lambdaLLR*llrObj + params.lambdaTemporalFT*FTObjT;

end

function g = grad(x,y,E,Et,params,shift_idx)%***********************************************

l1Smooth = 1e-15;
[ny,nx,K] = size(x);
TVOPY = TV_Temp_AnyDim(1);
TVOPX = TV_Temp_AnyDim(2);
TVOPT = TV_Temp_AnyDim(3);
% FTOPT = FT_Temp_AnyDim(3);
ny_diad = 2^nextpow2(ny);

% L2-norm part
L2Grad = 2.*(Et(E(x)-y));
% L2Grad = (Et(E(x)-y));

WavGrad=0;
if params.lambdaWav > 0
    WavGrad = zeros(size(x));
    parfor meas = 1:K
        WavOP=Waveletxyt('Daubechies',4,6,[ny,nx]);
        w = WavOP*zpad(x(:,:,meas),[ny_diad ny_diad]);
        WavGrad(:,:,meas)=crop(WavOP'*(w.*(w.*conj(w)+l1Smooth).^(-0.5)),[ny nx]);
        %         w = WavOP*x(:,:,meas);
        %         WavGrad(:,:,meas) = WavOP'*(w.*(w.*conj(w)+l1Smooth).^(-0.5));
    end
end

TVGradX = 0;
if params.lambdaSpatialTV > 0
    w = TVOPX*x;
    TVGradX = reshape(TVOPX'*(w.*(w.*conj(w)+l1Smooth).^(-0.5)),[ny nx K]);
end

TVGradY = 0;
if params.lambdaSpatialTV > 0
    w = TVOPY*x;
    TVGradY = reshape(TVOPY'*(w.*(w.*conj(w)+l1Smooth).^(-0.5)),[ny nx K]);
end

TVGradT = 0;
if params.lambdaTemporalTV > 0
    if ~isempty(params.Phi)
        Nex=size(params.Phi,1);
        I = reshape(reshape(x,[ny*nx,K])*params.Phi',[ny,nx,Nex]);
        
        %     w = TVOPT*x;
        w = TVOPT * I;
        %     TVGradT = reshape(TVOPT'*(w.*(w.*conj(w)+l1Smooth).^(-0.5)),[ny nx K]);
        TVGradT = reshape(TVOPT'*(w.*(w.*conj(w)+l1Smooth).^(-0.5)),[ny nx Nex]);
        TVGradT = reshape(reshape(TVGradT,[ny*nx,Nex])*params.Phi,[ny,nx,K]);
    else
        w = TVOPT*x;
        TVGradT = TVOPT'*(w.*(w.*conj(w)+l1Smooth).^(-0.5));
    end
end

FTGradT = 0;
if params.lambdaTemporalFT > 0
    w = FTOPT*x;
    FTGradT = reshape(FTOPT'*(w.*(w.*conj(w)+l1Smooth).^(-0.5)),[ny nx K]);
end

llrGrad = 0;
if params.lambdaLLR > 0
    NN = ceil(ny / params.block_dim(1)) * params.block_dim(1);
    [w,U_big,V_big,shift_idx] = imageToLowRankPatches(zpad(x,NN,NN,K), params.block_dim, params.block_step,shift_idx);
    w = (w.*(w.*conj(w)+l1Smooth).^(-0.5));
    [ llrGrad ] = crop( lowRankPatchesToImage([NN NN K], params.block_dim, params.block_step, shift_idx, U_big,V_big,w),ny,nx,K);
end

g = L2Grad + params.lambdaSpatialTV*TVGradX + params.lambdaSpatialTV*TVGradY + ...
    params.lambdaWav*WavGrad + params.lambdaTemporalTV*TVGradT + params.lambdaLLR*llrGrad + ...
    params.lambdaTemporalFT*FTGradT;
end

