function [out] = WASABI_vis(df0_lims, db1_lims, mode, ImagesNorm, df0_Map, db1_Map, C_Map, D_Map, R2_Map, mask_fit, B1, f_off, tau)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% create colormaps
mymap1 = get_cmp('blue_red', 1000);
mymap1(1,:) = 0;
df0_Map(mask_fit==0) = -1e6;
db1_Map(mask_fit==0) = -1e6;

mymap2 = jet(1024);
mymap2(1,:) = 0;

%% vis
if mode==1
    figure()
    ax(1) = subplot(1,3,1);
    imagesc(df0_Map, df0_lims);  axis image; axis off; colormap(ax(1), mymap1); colorbar; title('B0 map')
    ax(2) = subplot(1,3,2);
    imagesc(db1_Map, db1_lims); axis image; axis off; colormap(ax(2), mymap1); colorbar; title('B1+ map')
    ax(3) = subplot(1,3,3);
    imagesc(R2_Map,[0.9,1]);    axis image; axis off; colormap(ax(3), mymap2); colorbar; title('R2 map')
end

%% interactive visulaization
if mode==2
    gamma  = 2.67522/2/pi * 1e8;
    xdata  = linspace(f_off(1), f_off(end), 1000);
    FITfun = @(P) abs(P(1)-P(2)*(sin(atan((gamma.*P(3))./(xdata-P(4))))).^2 .* (sin((((gamma.*P(3)).^2+(xdata-P(4)).^2).^(1/2)).*(2*pi*tau/2))).^2);
     
    figure()
    ax(1) = subplot(2,2,1);
    imagesc(df0_Map, df0_lims);  axis image; axis off; colormap(ax(1), mymap1); colorbar; title('B0 map')
    ax(2) = subplot(2,2,2);
    imagesc(db1_Map, db1_lims); axis image; axis off; colormap(ax(2), mymap1); colorbar; title('B1+ map')
    
    vis = 1;
    while vis==1
        subplot(2,2,1)
        [x,y] = ginput(1);
        if x<10
            vis=0;
        end
        x = round(x);
        y = round(y);    
        temp_ydata     = ImagesNorm(:,y,x);
        temp_P         = [C_Map(y,x) D_Map(y,x) db1_Map(y,x)*B1 df0_Map(y,x)];
        temp_ydata_fit = FITfun(temp_P);    
        ax(3) = subplot(2,2,[3 4]);
        cla()
        hold on        
        plot(xdata, temp_ydata_fit, 'k-', 'LineWidth', 2)
        plot(f_off', temp_ydata, 'r.', 'MarkerSize', 20)
        sgtitle(['df0: ' num2str(df0_Map(y,x), '%.1f') 'Hz   dB1: ' num2str(db1_Map(y,x), '%.2f')]) 
        hold off
    end
end

out = [];

end

