function out = MRF_cardio_check_invivo_timings(SEQ, study_path, twix_obj, PULSEQ)

    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

    % number of noise pre-scans
    if isfield(PULSEQ.SPI, 'Nnoise')
        Nnoise = PULSEQ.SPI.Nnoise;
    else
        Nnoise = 0;
    end

    % read ECG monitoring from PMU data
    ecg = read_PMU_ECG(study_path);
    
    % time axis of PMU: 2.5ms raster
    tPMU = double(ecg.timestamp') * 0.0025;
    tPMU_adc_first = double(twix_obj.image.timestamp(1+Nnoise)) * 0.0025; % PMU time of first adc
    
    % time axis of SEQ: 1us raster
    tSEQ = (1:numel(SEQ.FULL.RF))' * SEQ.dt;
    tSEQ_adc_first = tSEQ(min(find(SEQ.FULL.ADC_ON_OFF==1))); % SEQ time of first adc
    
    % shift time of PMU axis
    tECG = tPMU - tPMU_adc_first + tSEQ_adc_first;
    clear tPMU tPMU_adc_first tSEQ_adc_first;
    
    % ECG signals
    ECG = double(ecg.filteredData) / max(max(double(ecg.filteredData)));
    clear ecg;
    
    % plot ECG and sequence
    fig = figure();
    ax1 = subplot(3,1,1);
    hold on
    plot(tECG, ECG(1,:), '-', 'color', 0.0*[1 1 1], 'LineWidth', 1);
    plot(tECG, ECG(2,:), '-', 'color', 0.1*[1 1 1], 'LineWidth', 0.1);
    plot(tECG, ECG(3,:), '-', 'color', 0.2*[1 1 1], 'LineWidth', 0.1);
    plot(tECG, ECG(4,:), '-', 'color', 0.3*[1 1 1], 'LineWidth', 0.1);
    for j=1:numel(SEQ.trig_timings)
        xline(SEQ.trig_timings(j), 'g-', 'LineWidth', 2)
    end
    xlim([tECG(1) tECG(end)]);
    ylabel('ecg [a.u.]');
    yticks([1]);
    set(gca, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
    
    ax2 = subplot(3,1,2);
    plot(tSEQ(SEQ.FULL.ADC_ON_OFF==1), SEQ.FULL.ADC_ON_OFF(SEQ.FULL.ADC_ON_OFF==1), 'rx');
    xlim([tSEQ(1) tSEQ(end)]);
    ylabel('ADC on/off');
    yticks([1]);
    set(gca, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
    
    ax3 = subplot(3,1,3);
    plot(tSEQ, abs(SEQ.FULL.RF), 'b-');
    xlim([tSEQ(1) tSEQ(end)]);
    ylabel('RF magnitude [Hz]');
    set(gca, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
    
    linkaxes([ax1 ax2 ax3], 'x')
    h = zoom(fig);
    setAxesZoomMotion(h, ax1, 'horizontal');
    clear ax1 ax2 ax3 h fig;
    
    out = [];
    
end