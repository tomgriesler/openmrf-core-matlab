function [FAs, TRs] = MRF_get_FAs_TRs(pattern_name, plotflag)

    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

    pattern_path = pulseq_get_path('MRF_find_patterns');
    if isempty(pattern_name)
    
    else
        pattern_path = [pattern_path pattern_name];
    end
    
    pattern_path_fa = [pattern_path '/fa.txt'];
    pattern_path_tr = [pattern_path '/tr.txt'];
    temp_ID         = fopen(pattern_path_fa,'r');
    FAs             = fscanf(temp_ID, '%f');
    fclose(temp_ID);
    temp_ID         = fopen(pattern_path_tr,'r');
    TRs             = fscanf(temp_ID, '%f');
    fclose(temp_ID);
    
    FAs = FAs *pi/180; % [rad]
    TRs = TRs *1e-3;   % [s]
    
    if plotflag==1
        figure()
        subplot(1,2,1)
        plot(FAs*180/pi, '.-')
        xlabel('number of acquisitions')
        ylabel('flip angles [deg]')
        axis square
        set(gca, 'FontName', 'Arial')
        subplot(1,2,2)
        plot(TRs*1e3, '.-')
        xlabel('number of acquisitions')
        ylabel('repetition time [ms]')
        axis square
        set(gca, 'FontName', 'Arial')
        sgtitle(['import MRF pattern: ' pattern_name])
    end

end