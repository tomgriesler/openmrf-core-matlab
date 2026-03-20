function ADIASL = ADIASL_init(ADIASL, FOV, system)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% adiabatic spin-lock preparation: global inversion pulse
if strcmp(ADIASL.mode, 'on')

    % set default parameters
    if ~isfield(ADIASL, 'phi')
        ADIASL.phi = 0;
    end
    if ~isfield(ADIASL, 'tau')
        ADIASL.tau = 0.015;
    end
    if ~isfield(ADIASL, 'f1max')
        ADIASL.f1max = 500;
    end
    if ~isfield(ADIASL, 'beta')
        ADIASL.beta = 5;
    end
    if ~isfield(ADIASL, 'dwmax')
        ADIASL.dfmax = 200;
    end
    if ~isfield(ADIASL, 'phi_list')
        ADIASL.phi_list = repmat([0 1 1 0   1 0 0 1   0 0 1 1   1 1 0 0], 1, 100)';
    end

    % create dummy pulse objects
    ADIASL.rf = mr.makeSincPulse( pi/2, system, 'duration', ADIASL.tau, 'use', 'other');
    ADIASL.rf.signal = ADIASL.rf.signal * 0;

    % calculate waveform of hyperbolic sechans:
    ADIASL.f1       = ADIASL.f1max * sech(ADIASL.beta*(2*ADIASL.rf.t/ADIASL.tau-1));
    ADIASL.df_down  = -2*ADIASL.dfmax*tanh(ADIASL.beta*(2*ADIASL.rf.t/ADIASL.tau-1)); % !!!
    ADIASL.df_up    = -2*ADIASL.dfmax*tanh(ADIASL.beta*(2*ADIASL.rf.t/ADIASL.tau-1)); % improved performance if both pulses are identical; compare to Coletti et al!
    ADIASL.phi_down = cumsum(2*pi*ADIASL.df_down) * system.rfRasterTime;
    ADIASL.phi_up   = cumsum(2*pi*ADIASL.df_up)   * system.rfRasterTime;
    ADIASL.phi_down = wrapTo2Pi(ADIASL.phi_down - ADIASL.phi_down(round(numel(ADIASL.f1)/2)));
    ADIASL.phi_up   = wrapTo2Pi(ADIASL.phi_up   - ADIASL.phi_up(  round(numel(ADIASL.f1)/2)));

    % calc prep times
    ADIASL.adia_prep_times = ADIASL.N_HS * mr.calcDuration(ADIASL.rf);

    % calculate crusher objects
    if ~isfield(ADIASL, 'crush_nTwists_x')
        ADIASL.crush_nTwists_x = 4;   % [] number of 2pi twists in x direction
    end
    if ~isfield(ADIASL, 'crush_nTwists_y')
        ADIASL.crush_nTwists_y = 4;   % [] number of 2pi twists in y direction
    end
    if ~isfield(ADIASL, 'crush_nTwists_z')
        ADIASL.crush_nTwists_z = 11.3;   % [] number of 2pi twists in z direction
    end
    if ~isfield(ADIASL, 'crush_lim_grad')
        ADIASL.crush_lim_grad = 1/sqrt(3); % reduce crusher gradient amplitude from nominal limit
    end    
    if ~isfield(ADIASL, 'crush_lim_slew')
        ADIASL.crush_lim_slew = 1/sqrt(3); % reduce crusher gradient slew rate from nominal limit to avoid stimulation
    end
    [ADIASL.gx_crush, ADIASL.gy_crush, ADIASL.gz_crush] = CRUSH_x_y_z(ADIASL.crush_nTwists_x, ADIASL.crush_nTwists_y, ADIASL.crush_nTwists_z, FOV.dx, FOV.dy, FOV.dz, ADIASL.crush_lim_grad, ADIASL.crush_lim_slew, system);
    ADIASL.tcrush = mr.calcDuration(ADIASL.gx_crush, ADIASL.gy_crush, ADIASL.gz_crush);
    
end

end