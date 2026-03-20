function [db1_Map, df0_Map, R2_Map, C_Map, D_Map] = WASABI_fit(ImageArr, mask_fit, f_off, tau, b1_start)   

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% ----- input -----
% ImageArr:   magnitude images
% mask_fit:   binary fit mask
% f_off:      [Hz] frequency offsetes of wasabi pulses
% tp:         [s]  duration of wasabi pulses
% b1_start:   [T] nominal amplitude of wasabi pulses

% ----- output -----
% dB1_Map % [0.8 ... 1.2] b1+ map
% df0_Map % [Hz] water shift map
% R2_Map  % [0...1] coefficient of determination
% C_Map   % output for fit control
% D_Map   % output for fit control

%% WASABI fit:
% Schuenke P et al. Magn Reson Med. 2017 Feb;77(2):571-580. doi: 10.1002/mrm.26133
% fit equation:  Z(dw) = abs(c-d*(sin(atan((gamma.*B1)./(x-dw)))).^2 .* (sin((((gamma.*B1).^2+(x-dw).^2).^(1/2)).*(2*pi*tp/2))).^2)
gamma  = 2.67522/2/pi * 1e8;
FITfun = @(P) abs(P(1)-P(2)*(sin(atan((gamma.*P(3))./(f_off-P(4))))).^2 .* (sin((((gamma.*P(3)).^2+(f_off-P(4)).^2).^(1/2)).*(2*pi*tau/2))).^2);

%% init
Ny = size(ImageArr,2); % Nxy
Nx = size(ImageArr,3); % Nxy

df0_Map = zeros(Nx,Nx); % [Hz]
db1_Map = zeros(Ny,Nx); % [0.8 ... 1.2]
C_Map   = zeros(Ny,Nx);
D_Map   = zeros(Ny,Nx);
R2_Map  = zeros(Ny,Nx);

%% create fit mask
if isempty(mask_fit)
    mask_fit = get_mask_fit(squeeze(mean(abs(ImageArr))));
end

%% prefit with lookup table
look_up_c    = [0.2 : 0.1 : 1]';
look_up_d    = [0.5 : 0.25 : 2.0]';
look_up_b1   = [0.65 : 0.05 : 1.35]' * b1_start;
look_up_df0  = [-50 : 10 : -10, -7.5 : 2.5 : 7.5, 10 : 10 : 50]';
[look_up_c, look_up_d, look_up_b1, look_up_df0] = ndgrid(look_up_c, look_up_d, look_up_b1, look_up_df0);
look_up_c    = look_up_c(:);
look_up_d    = look_up_d(:);
look_up_b1   = look_up_b1(:);
look_up_df0  = look_up_df0(:);
n_dict       = numel(look_up_df0);
look_up_data = zeros(n_dict,numel(f_off));

for j = 1:n_dict
    P = [look_up_c(j), look_up_d(j), look_up_b1(j), look_up_df0(j)];
    temp = FITfun(P);
    look_up_data(j,:) = temp;
end

%% wasabi mapping
parfor y=1:Ny
for    x=1:Nx
    if mask_fit(y,x) == 1

    ydata        = squeeze(ImageArr(:,y,x));
    look_up_diff = sum((look_up_data' - repmat(ydata, 1, n_dict)).^2);
    look_up_ind  = find(look_up_diff==min(look_up_diff));
    look_up_ind  = look_up_ind(1);
    P_start      = [look_up_c(look_up_ind) look_up_d(look_up_ind) look_up_b1(look_up_ind) look_up_df0(look_up_ind)];
    
    % fminsearch()
    options      = optimset('Display','off', 'TolFun', 1e-9, 'TolX', 1e-9);
    RSSfun       = @(P) sum( (ydata' - FITfun(P)).^2 ); 
    P_fit        = fminsearch(RSSfun, P_start, options);
    TSS          = sum( (ydata - mean(ydata)).^2 );
    R2_Map(y,x)  = 1 - RSSfun(P_fit) / TSS;
    C_Map(y,x)   = P_fit(1);
    D_Map(y,x)   = P_fit(2);
    db1_Map(y,x) = P_fit(3) / b1_start;
    df0_Map(y,x) = P_fit(4);

    end
end
end



end

