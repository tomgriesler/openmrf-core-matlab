function [ T12p_Map, S0_Map, R2_Map, RSS_Map ] = mg_map_T12p( ImageArr, tSL, mask_fit, mod_fit )

% Version: Maximilian Gram, 21.03.2024

% fit model:
% S(tSL) = S0 * exp(-tSL/T1p)
% ImageArr = abs(ImageArr);

%% default parameters
[NtSL, Ny, Nx] = size(ImageArr);

if isempty(mask_fit)
    mask_fit = ones(Ny, Nx);
end
if nargin<4
    mod_fit = [];
end
if isempty(mod_fit)
    mod_fit = 1; %  0: toolbox   1: fminsearch   2: fmincon
end

%% save results in arrays
T12p_Map = zeros(Ny, Nx);
S0_Map   = zeros(Ny, Nx);
R2_Map   = zeros(Ny, Nx);
RSS_Map  = zeros(Ny, Nx);

%% ----- prefit: find start values for S0 and T1p -----
ydata_pre = zeros(NtSL,1);
for j=1:NtSL
    temp1(:,:)     = ImageArr(j,:,:);
    temp2(:)       = temp1(mask_fit==1);
    ydata_pre(j,1) = mean(temp2);
end
clear temp1 temp2;

temp1(:,:)            = ImageArr(1,:,:);
temp2(:)              = temp1(mask_fit==1);
S0_start              = mean(temp2);
T1p_start             = tSL(end)*0.75;
[S0_start, T1p_start] = mg_fit_exp(tSL, ydata_pre, S0_start, T1p_start, mod_fit);
clear temp1 temp2
        
%% pixel:pixel fit (parallel)
parfor y=1:Ny
for x=1:Nx
    if (mask_fit(y,x) == 1)
        xdata               = tSL;
        ydata               = ImageArr(:,y,x);
        [S0, T12p, R2, RSS] = mg_fit_exp(xdata, ydata, S0_start, T1p_start, mod_fit);
        S0_Map(y,x)         = S0;
        T12p_Map(y,x)       = T12p;
        R2_Map(y,x)         = R2;
        RSS_Map(y,x)        = RSS;
    end
end
end

end

