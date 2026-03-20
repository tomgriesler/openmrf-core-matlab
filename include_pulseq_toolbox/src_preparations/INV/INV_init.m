function INV = INV_init(INV, FOV, system)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% magnetization inversion: global inversion pulse

    if ~isfield(INV,'rf_type')
        INV.rf_type = 'HYPSEC_inversion'; % default
    end    
   
    if strcmp(INV.rf_type, 'HYPSEC_inversion')
        if ~isfield(INV, 'tExc')
            INV.tExc = 10 *1e-3;  % [s] pulse duration
        end
        if ~isfield(INV, 'mu')
            INV.mu = 4.9;         % [ ] determines amplitude of frequency sweep
        end
        if ~isfield(INV, 'beta')
            INV.beta = 600;       % [Hz] peak amplitude
        end        
        [INV.rf] = SIGPY_HYPSEC(INV.beta, INV.mu, INV.tExc, 0, system);

    elseif strcmp(INV.rf_type, 'sinc')
        INV.tExc = 6 *1e-3;  % [s] pulse duration
        [INV.rf] = mr.makeSincPulse( pi, system, 'Duration', INV.tExc, 'SliceThickness', 1e6, 'use', 'inversion');

    elseif strcmp(INV.rf_type, 'sigpy_SLR')
        INV.tExc = 6 *1e-3;  % [s] pulse duration
        [INV.rf] = SIGPY_SLR(pi, INV.tExc, 0, 4, 'inv', 'ls', 0.01, 0.01, 0 ,1e6, system);

    elseif strcmp(INV.rf_type, 'AFP_inversion')
        ahp.tau      = 3 *1e-3;     % [s] half pulse duration
        ahp.wmax     = 600 * 2*pi;  % [rad/s]
        ahp          = AHP_get_params( ahp.tau, ahp.wmax);
        ahp.rf       = AFP_pulse( system, ahp.tau, ahp.wmax, ahp.ratio, ahp.beta1, ahp.beta2, 0 );
        INV.rf       = ahp.rf;
        INV.tExc     = 2*ahp.tau;
        clear ahp;
    end

    % calculate crusher objects
    if ~isfield(INV, 'crush_nTwists_x')
        INV.crush_nTwists_x = 4;   % [] number of 2pi twists in x direction
    end
    if ~isfield(INV, 'crush_nTwists_y')
        INV.crush_nTwists_y = 4;   % [] number of 2pi twists in y direction
    end
    if ~isfield(INV, 'crush_nTwists_z')
        INV.crush_nTwists_z = 11.6;   % [] number of 2pi twists in z direction
    end
    if ~isfield(INV, 'crush_lim_grad')
        INV.crush_lim_grad = 1/sqrt(3); % reduce crusher gradient amplitude from nominal limit
    end    
    if ~isfield(INV, 'crush_lim_slew')
        INV.crush_lim_slew = 1/sqrt(3); % reduce crusher gradient slew rate from nominal limit to avoid stimulation
    end
    [INV.gx_crush, INV.gy_crush, INV.gz_crush] = CRUSH_x_y_z(INV.crush_nTwists_x, INV.crush_nTwists_y, INV.crush_nTwists_z, FOV.dx, FOV.dy, FOV.dz, INV.crush_lim_grad, INV.crush_lim_slew, system);
    INV.tcrush = mr.calcDuration(INV.gx_crush, INV.gy_crush, INV.gz_crush);

    % recovery delay
    if isfield(INV, 'inv_rec_time')
        if size(INV.inv_rec_time,2) > size(INV.inv_rec_time,1)
            INV.inv_rec_time = INV.inv_rec_time';
        end
        for j=1:numel(INV.inv_rec_time)
            INV.inv_rec_delay(j,1) = mr.makeDelay( round(INV.inv_rec_time(j)/system.gradRasterTime)*system.gradRasterTime );
        end
    else
        INV.inv_rec_time  = system.gradRasterTime;
        INV.inv_rec_delay = mr.makeDelay( round(INV.inv_rec_time/system.gradRasterTime)*system.gradRasterTime );
    end

end