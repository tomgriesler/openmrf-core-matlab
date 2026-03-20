function [rawdata, noise, PULSEQ, study_info] = pulseq_read_meas(path_raw, path_backup, vendor)

    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
    
    % ----- Input: -----
    % path_raw:      file path to rawdata (.dat for Siemens or .mat for others)
    % path_backup:   file path to PULSEQ backup (.mat, not necessary for Siemens)
    % vendor:        'Siemens'  'GE'  'UnitedImaging'  'Philips'
    
    % ----- Output: -----
    % rawdata:     [coils x tr x adc] complex rawdata
    % noise:       [coils x rep x adc] complex noise pre-scans
    % PULSEQ       struct containing all sequence parameters
    % study_info:  struct containing scan/study specific header information

    switch vendor
        case 'Siemens' % for Siemens, we can automatically link rawdata and pulseq backups
            [twix_obj, study_info, PULSEQ] = pulseq_read_meas_siemens(path_raw);
            rawdata = SPI_get_rawdata(twix_obj);
            rawdata = permute(rawdata, [3, 1, 2]); % coils x tr x adc
            if ~isempty(path_backup)
                load(path_backup);
                warning('automatic PULSEQ backup was overwritten for Siemens scan!');
            end            
    
        case 'GE' % for all other vendors, we still need routines...
            load(path_raw);
            load(path_backup);
            study_info.time_stamps = time_stamps;
        
        case 'UnitedImaging'
            load(path_raw); 
            load(path_backup);
            study_info.time_stamps = time_stamps;
        
        case 'Philips'
            load(path_raw);
            load(path_backup);
            study_info.time_stamps = time_stamps;
    end

    if isfield(PULSEQ, 'SPI') && isfield(PULSEQ.SPI, 'Nnoise') && PULSEQ.SPI.Nnoise>0
        noise = rawdata(:,1:PULSEQ.SPI.Nnoise,:);
        rawdata(:,1:PULSEQ.SPI.Nnoise,:) = [];
    else
        noise = [];
    end

end
