function WASABI = WASABI_init(WASABI, FOV, system)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% calculate WASABI pulse object and S0 reference delay
WASABI.rf = mr.makeBlockPulse( pi/10, system, ...
                               'Duration', WASABI.tau, ...
                               'freqOffset', 0, ...
                               'PhaseOffset', WASABI.phase, ...
                               'use', 'preparation' );
WASABI.rf.signal = WASABI.rf.signal*0 + WASABI.f1;
WASABI.ref_delay = mr.calcDuration(WASABI.rf);

%% crusher gradients

% calculate crusher objects
if ~isfield(WASABI, 'crush_nTwists_x')
    WASABI.crush_nTwists_x = 4;   % [] number of 2pi twists in x direction
end
if ~isfield(WASABI, 'crush_nTwists_y')
    WASABI.crush_nTwists_y = 4;   % [] number of 2pi twists in y direction
end
if ~isfield(WASABI, 'crush_nTwists_z')
    WASABI.crush_nTwists_z = 11.3;   % [] number of 2pi twists in z direction
end
if ~isfield(WASABI, 'crush_lim_grad')
    WASABI.crush_lim_grad = 1/sqrt(3); % reduce crusher gradient amplitude from nominal limit
end    
if ~isfield(WASABI, 'crush_lim_slew')
    WASABI.crush_lim_slew = 1/sqrt(3); % reduce crusher gradient slew rate from nominal limit to avoid stimulation
end
[WASABI.gx_crush, WASABI.gy_crush, WASABI.gz_crush] = CRUSH_x_y_z(WASABI.crush_nTwists_x, WASABI.crush_nTwists_y, WASABI.crush_nTwists_z, FOV.dx, FOV.dy, FOV.dz, WASABI.crush_lim_grad, WASABI.crush_lim_slew, system);
WASABI.tcrush = mr.calcDuration(WASABI.gx_crush, WASABI.gy_crush, WASABI.gz_crush);

end

