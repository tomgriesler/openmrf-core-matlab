function [df0_Map, dB0_Map, R2_Map, PhaseCorr] = mg_map_df0(ImageArr, TEs, mask_fit, df0_lims, unwrap_flag, plot_flag)

% Version: Maximilian Gram, 21.03.2024

% input:
% ImageArr: complex MRI images
% TEs: echo times [s]
% mask_fit: binary mask
% df0_lims: limits [Hz] for colorbar
% unwrap_flag: [0 0] -> no unwrapping
%              [1 0] -> matlab unwrapping
%              [0 1] -> mgram unwrapping
%              [1 1] -> matlab & mgram unwrapping
% plot_flag:   0 -> nothing
%              1 -> show df0 and r2 maps
%              2 -> interactive tool

% output:
% df0_Map: offresonance map [Hz]
% dB0_Map: offresonance map [nT]
% R2_Map: fit quality
% PhaseCorr: corrected unwrapped phases

%% set defaults

if nargin<6 || isempty(plot_flag)
    plot_flag = 0; % no graphical output
end

if nargin<5 || isempty(unwrap_flag)
    unwrap_flag = [1 0]; % matlab unwrapper
end

if nargin<4
    df0_lims = [];
end

if nargin<3 || isempty(mask_fit)
    mask_fit = mg_get_mask_fit(squeeze(mean(abs(ImageArr))), 'holes');
end

if size(TEs,1) < size(TEs,2)
    TEs = TEs';
end

%% pixel:pixel fit
[nTE, ny, nx] = size(ImageArr);
df0_Map       = zeros(ny, nx);
R2_Map        = zeros(ny, nx);
offset_map    = zeros(ny, nx);
PhaseCorr     = zeros(nTE, ny, nx);

parfor y=1:ny
for    x=1:nx
    if mask_fit(y,x)==1

        P_corr           = mg_unwrap_phase(squeeze(ImageArr(:,y,x)), TEs, unwrap_flag);
        [A, B, R2]       = mg_fit_lin(TEs, P_corr);
        offset_map(y,x)  = A;
        df0_Map(y,x)     = B;
        R2_Map(y,x)      = R2;
        PhaseCorr(:,y,x) = P_corr(:);
        
    end
end
end

df0_Map = df0_Map/2/pi;         % [Hz]
dB0_Map = mg_f2B(df0_Map) *1e9; % [nT]

%% graphical output: 1
if plot_flag > 0
    if isempty(df0_lims)
        df0_lims = max(abs(df0_Map(mask_fit==1)));
        df0_lims = ceil(df0_lims/25)*25;
        df0_lims = [-df0_lims df0_lims];
    end
    df0_Map_ = df0_Map;
    df0_Map_(mask_fit==0) = -Inf;
    mycmp1 = get_cmp('blue_red', 1024);
    mycmp1(1,:) = 0;
    mycmp2 = jet(1024);
    mycmp2(1,:) = 0;
    figure()
    ax(1) = subplot(1,2,1);
    imagesc(df0_Map_, df0_lims); axis image; axis off; colormap(ax(1), mycmp1); colorbar; title('df0 Map [Hz]')
    ax(2) = subplot(1,2,2);
    imagesc(R2_Map, [0.9, 1]); axis image; axis off; colormap(ax(2), mycmp2); colorbar; title('R2 Map [0...1]')
end

%% graphical output: 2
if plot_flag > 1

    TEs_fit = linspace(TEs(1), TEs(end), 100);
    FITfun  = @(P) P(1) + TEs_fit * P(2) * 2*pi;
    
    figure()
    subplot(1,2,1)
    imagesc(df0_Map_, df0_lims); axis image; axis off; colormap(mycmp1); colorbar; title('df0 Map [Hz]')
    
    while 1==1    
        subplot(1,2,1)
        [x,y] = ginput(1);
        x = round(x);
        y = round(y);
        
        subplot(1,2,2)
        cla()
        hold on
        plot(TEs*1e3,     squeeze(PhaseCorr(:,y,x)), 'ro')
        plot(TEs*1e3,     squeeze(angle(ImageArr(:,y,x))), 'r.') 
        plot(TEs_fit*1e3, FITfun([offset_map(y,x), df0_Map_(y,x)]), 'b-')
        ylim([min([min(PhaseCorr(:)), -pi]), max([max(PhaseCorr(:)), pi])])
        xlim([TEs(1)*1e3 TEs(end)*1e3])
        xlabel(['echo time TE [ms]'])
        ylabel(['unwrapped phase [rad]'])
        title(['df0 = ' num2str(df0_Map_(y,x),'%.1f') 'Hz   R2 = ' num2str(R2_Map(y,x),'%.3f')])
        hold off
    end

end

end

