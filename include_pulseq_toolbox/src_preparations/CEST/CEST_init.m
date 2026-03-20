function CEST = CEST_init(CEST, FOV, system)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% calculate CEST pulse object and S0 reference delay
CEST.rf        = SL_get_SL_pulse(system, CEST.tau, CEST.f1, CEST.phase, 0);
CEST.ref_delay = mr.calcDuration(CEST.rf);

%% crusher gradients

if ~isfield(CEST, 'crush_nTwists_x')
    CEST.crush_nTwists_x = 4;   % [] number of 2pi twists in x direction
end
if ~isfield(CEST, 'crush_nTwists_y')
    CEST.crush_nTwists_y = 4;   % [] number of 2pi twists in y direction
end
if ~isfield(CEST, 'crush_nTwists_z')
    CEST.crush_nTwists_z = 11.3;   % [] number of 2pi twists in z direction
end
if ~isfield(CEST, 'crush_lim_grad')
    CEST.crush_lim_grad = 1/sqrt(3); % reduce crusher gradient amplitude from nominal limit
end    
if ~isfield(CEST, 'crush_lim_slew')
    CEST.crush_lim_slew = 1/sqrt(3); % reduce crusher gradient slew rate from nominal limit to avoid stimulation
end
[CEST.gx_crush, CEST.gy_crush, CEST.gz_crush] = CRUSH_x_y_z(CEST.crush_nTwists_x, CEST.crush_nTwists_y, CEST.crush_nTwists_z, FOV.dx, FOV.dy, FOV.dz, CEST.crush_lim_grad, CEST.crush_lim_slew, system);
CEST.tcrush = mr.calcDuration(CEST.gx_crush, CEST.gy_crush, CEST.gz_crush);

end

