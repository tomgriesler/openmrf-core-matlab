function [alpha_slice, alpha_med, alpha_mean, alpha_std, corr_scale] = mg_slice_sim(rf, dz, gz, alpha, dt, N_iso, sfac)

    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
    
    %% bloch simulation for isochromats
    z     = linspace(-1/2, 1/2, N_iso) *dz *sfac;
    df    = z * gz;
    f1    = abs(rf);
    phi   = angle(rf);
    M_iso = zeros(N_iso, 3);
    
    parfor j=1:N_iso    
        dw_ = 2*pi*df(j);
        M_  = [0; 0; 1];    
        for k=1:numel(rf)
            w1x_ = 2*pi * f1(k) * cos(phi(k));
            w1y_ = 2*pi * f1(k) * sin(phi(k));
            B_   = [  0     dw_    -w1y_;
                     -dw_   0      w1x_;
                     w1y_   -w1x_  0 ];
            M_ = expm(B_*dt) * M_;
        end    
        M_iso(j,:) = M_(:);
    end
    
    %% calculate slice profile of isochromats
    Mz          = M_iso(:,3);
    Mxy         = sqrt( M_iso(:,1).^2 + M_iso(:,2).^2 );
    alpha_iso   = real(heaviside(sign(Mz)) .* asin(Mxy) + heaviside(-sign(Mz)) .* (pi-asin(Mxy)));
    
    %% show results
    temp = z*0;
    temp(abs(z)<dz/2) = 1;
    temp = alpha_iso(temp==1);
    alpha_med  = median(temp);
    alpha_mean = mean(temp);
    alpha_std  = std(temp);
    corr_scale = alpha / alpha_med;
    
    figure()
    hold on
    xline(-dz/2*1e3, 'k--', 'LineWidth', 2)
    xline( dz/2*1e3, 'k--', 'LineWidth', 2)
    yline(alpha*180/pi,     'g--', 'LineWidth', 2, 'Label', 'nominal')
    yline(alpha_med*180/pi, 'r--', 'LineWidth', 2, 'Label', 'median')
    plot( z*1e3, alpha_iso*180/pi, 'b-', 'LineWidth', 3 )
    errorbar(-0.9*max(z)*1e3, alpha_mean*180/pi, alpha_std*180/pi, 'ko', 'LineWidth', 2)
    title(['nom: '     num2str(alpha*180/pi,'%.1f'), ...
           '  med: '   num2str(alpha_med*180/pi,'%.1f'), ...
           '  mean: '  num2str(alpha_mean*180/pi,'%.1f') '+-' num2str(alpha_std*180/pi,'%.1f'), ...
           '  corr: '  num2str(corr_scale*100-100,'%.1f') '%'])
    xlim([-1 1]*max(z*1e3))
    xlabel('slice profile z [mm]')
    ylabel('flip angle [deg]')
    set(gca,'linewidth', 2, 'fontsize', 12, 'fontname', 'arial', 'fontweight', 'bold')
    hold off

end

