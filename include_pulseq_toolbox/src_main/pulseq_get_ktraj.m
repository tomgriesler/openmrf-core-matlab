function [ktraj_adc, ktraj_full] = pulseq_get_ktraj(seq, flag_plot)
% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
    warning('OFF', 'mr:restoreShape');
    [ktraj_adc, ~, ktraj_full] = seq.calculateKspacePP();
    if flag_plot==1
        figure();
        hold on
        plot( ktraj_full(1,:), ktraj_full(2,:), 'b-');
        plot( ktraj_adc(1,:),  ktraj_adc(2,:), 'r.');
        axis('equal');
        xlabel('kx [1/m]');
        ylabel('ky [1/m]');
        title('k-space trajectory');
    elseif flag_plot==2
        figure();
        hold on
        kmax = max(abs(ktraj_full(:)));
        plot3([-1 1]*kmax, [0 0], [0 0], 'k-')
        plot3([0 0], [-1 1]*kmax, [0 0], 'k-')
        plot3([0 0], [0 0], [-1 1]*kmax, 'k-')
        plot3( ktraj_full(1,:), ktraj_full(2,:), ktraj_full(3,:), 'b-');
        plot3( ktraj_adc(1,:),  ktraj_adc(2,:), ktraj_adc(3,:), 'r.');
        axis('equal');
        xlabel('kx [1/m]');
        ylabel('ky [1/m]');
        zlabel('kz [1/m]');
        title('k-space trajectory');
        view([-37.5 -30]); 
   end

end