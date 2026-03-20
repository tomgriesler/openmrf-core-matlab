function [BSSFP, ktraj_adc, ktraj_full] = BSSFP_init(BSSFP, FOV, system)

% ---------------------------------------------------------
% ---------- init parameters and pulseq objects -----------
% ------------------- readout: 3D BSSFP -------------------
% ---------------------------------------------------------

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% limit max grad strength and slew rate
if ~isfield(BSSFP, 'max_grad')
    BSSFP.max_grad = 1.0;
end
if ~isfield(BSSFP, 'max_slew')
    BSSFP.max_slew = 1.0;
end
system.maxGrad = system.maxGrad * BSSFP.max_grad;
system.maxSlew = system.maxSlew * BSSFP.max_slew;

%% calculate rf objects
BSSFP.n_FAs = numel(BSSFP.exc_flipangle);  % allow different FAs for ramp up
for j = 1 : BSSFP.n_FAs
    if strcmp(BSSFP.exc_mode, 'sinc')
        [BSSFP.rf(j), BSSFP.gz] = mr.makeSincPulse( BSSFP.exc_flipangle(j), ...
                                                    system, ...
                                                    'Duration', BSSFP.exc_time, ...
                                                    'timeBwProduct', BSSFP.exc_tbw, ...
                                                    'apodization', 0.5, ...
                                                    'SliceThickness', FOV.dz, ...
                                                    'use', 'excitation'); 
    elseif strcmp(BSSFP.exc_mode, 'sigpy_SLR')
        [BSSFP.rf(j), BSSFP.gz] = SIGPY_SLR( BSSFP.exc_flipangle(j), BSSFP.exc_time, 0, BSSFP.exc_tbw, BSSFP.exc_shape, 'ls', 0.01, 0.01, 0 , FOV.dz, system);
    end
    BSSFP.rf(j).freqOffset  = BSSFP.gz.amplitude * FOV.z_offset;
    BSSFP.rf(j).phaseOffset = 0; % changed by phase cycling
end
clear j;

%% k-space parameters
BSSFP.delta_kx = 1/FOV.fov_x;
BSSFP.delta_ky = 1/FOV.fov_y;
BSSFP.delta_kz = 1/FOV.fov_z;
BSSFP.kx_area  = FOV.Nx * BSSFP.delta_kx;
BSSFP.ky_area  = (-FOV.Ny/2 : 1 : FOV.Ny/2-1)' * BSSFP.delta_ky;
BSSFP.kz_area  = (-FOV.Nz/2 : 1 : FOV.Nz/2-1)' * BSSFP.delta_kz;

%% calculate adc and read gradient
if BSSFP.os_mode==0
    BSSFP.adcNSamples = FOV.Nx;
elseif BSSFP.os_mode==1
    BSSFP.adcNSamples = 2*FOV.Nx;
end
BSSFP.adcNSamples = ceil(BSSFP.adcNSamples/64)*64;
BSSFP.adcDwell    = round( BSSFP.t_acq / BSSFP.adcNSamples / system.adcRasterTime ) * system.adcRasterTime;
BSSFP.adc         = mr.makeAdc(BSSFP.adcNSamples, 'Dwell', BSSFP.adcDwell); 
BSSFP.adcBW       = 1 / BSSFP.adcDwell;
BSSFP.adcDur      = mr.calcDuration(BSSFP.adc);
BSSFP.t_acq       = ceil( BSSFP.adcDur / system.blockDurationRaster ) * system.blockDurationRaster;
BSSFP.gx          = mr.makeTrapezoid('x', 'FlatTime', BSSFP.t_acq, 'amplitude', BSSFP.kx_area/BSSFP.adcDur, 'system', system);
BSSFP.adc.delay   = BSSFP.gx.riseTime + (BSSFP.t_acq - BSSFP.adcDur) / 2;
BSSFP.adc.delay   = round( BSSFP.adc.delay / system.adcRasterTime ) * system.adcRasterTime;

%% calculate prephasers & rewinders

if ~isfield(BSSFP, 't_pre_rew')
    error('define BSSFP.t_pre_rew: this is the duration for xyz prephaser and rewinder gradients!');
end
BSSFP.t_pre_rew = round(BSSFP.t_pre_rew/system.blockDurationRaster) * system.blockDurationRaster;

% x prephaser & rewinder
BSSFP.gx_pre_rew = mr.makeTrapezoid('x', 'Area', -BSSFP.gx.area/2, 'Duration', BSSFP.t_pre_rew, 'system', system);

% y prephaser & rewinder
for j = 1:FOV.Ny
    BSSFP.gy_pre(j) = mr.makeTrapezoid('y', 'Area',  BSSFP.ky_area(j), 'Duration', BSSFP.t_pre_rew, 'system', system);
    BSSFP.gy_rew(j) = mr.makeTrapezoid('y', 'Area', -BSSFP.ky_area(j), 'Duration', BSSFP.t_pre_rew, 'system', system);
end

% y prephaser & rewinder
for j = 1:FOV.Nz
    BSSFP.gz_pre(j) = mr.makeTrapezoid('z', 'Area',  BSSFP.kz_area(j) - BSSFP.gz.area/2, 'Duration', BSSFP.t_pre_rew, 'system', system);
    BSSFP.gz_rew(j) = mr.makeTrapezoid('z', 'Area', -BSSFP.kz_area(j) - BSSFP.gz.area/2, 'Duration', BSSFP.t_pre_rew, 'system', system);
end

%% calculate filling delays for TR: we fix TE = TR/2
tr_min = mr.calcDuration(BSSFP.rf, BSSFP.gz) + ...
         mr.calcDuration(BSSFP.gx_pre_rew, BSSFP.gy_pre(1), BSSFP.gz_pre(1)) + ...
         mr.calcDuration(BSSFP.gx, BSSFP.adc) + ...
         mr.calcDuration(BSSFP.gx_pre_rew, BSSFP.gy_rew(1), BSSFP.gz_rew(1));

if ~isfield(BSSFP, 'TR')
    BSSFP.TR = [];
end
if isempty(BSSFP.TR)
    BSSFP.TR = 0;
end

BSSFP.tr_fill = (BSSFP.TR - tr_min) / 2;
BSSFP.tr_fill = round(BSSFP.tr_fill / system.blockDurationRaster);
if BSSFP.tr_fill < 1
    BSSFP.tr_fill = 1;
end
BSSFP.tr_fill = mr.makeDelay(BSSFP.tr_fill * system.blockDurationRaster);
BSSFP.TR      = tr_min + mr.calcDuration(BSSFP.tr_fill)*2;
BSSFP.TE      = BSSFP.TR / 2;

%% calculate kspace trajectory
ktraj_adc  = [];
ktraj_full = [];

% only use for debugging
% seq = mr.Sequence(system);
% for loop_z = 1:FOV.Nz
%     for loop_y = 1:FOV.Ny
%         BSSFP_add();
%     end
% end
% [ktraj_adc, ktraj_full] = pulseq_get_ktraj(seq, 1);

end
