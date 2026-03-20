function [TRAJ, FOV] = TRAJ_init(TRAJ, SPI, system)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% set default values
if ~isfield(TRAJ, 'method')
    TRAJ.method = 'duyn';
end
if ~isfield(TRAJ, 'Nav')
    TRAJ.Nav = 10;
end
if ~isfield(TRAJ, 'Ndummy')
    TRAJ.Ndummy = 50;
end
if ~isfield(TRAJ, 'Trec')
    TRAJ.Trec = 150 *1e-3;
end
if ~isfield(TRAJ, 'slice_thickness')
    TRAJ.slice_thickness = 2 *1e-3;
end
if ~isfield(TRAJ, 'slice_offset')
    TRAJ.slice_offset = 50 *1e-3;
end
if ~isfield(TRAJ, 'exc_time')
    TRAJ.exc_time = 4 *1e-3;
end
if ~isfield(TRAJ, 'exc_tbw')
    TRAJ.exc_tbw = 4;
end
if ~isfield(TRAJ, 'exc_fa')
    TRAJ.exc_fa = 20 *pi/180;
end
if ~isfield(TRAJ, 'spoil_nTwists')
    TRAJ.spoil_nTwists = 8;
end
if ~isfield(TRAJ, 'lim_grad')
    TRAJ.lim_grad = 0.75;
end
if ~isfield(TRAJ, 'lim_slew')
    TRAJ.lim_slew = 0.75;
end

%% create adc for trajectory measurement and rotate projections
TRAJ.NR  = numel(SPI.rot);
TRAJ.adc = SPI.adc;

temp_g = [SPI.gx.waveform, SPI.gy.waveform, SPI.gz.waveform];
for j=1:TRAJ.NR
    temp_g_rot = (mr.aux.quat.toRotMat(SPI.rot(j).rotQuaternion) * temp_g')';
    TRAJ.gx(j,1) = SPI.gx;
    TRAJ.gy(j,1) = SPI.gy;
    TRAJ.gx(j,1).waveform = temp_g_rot(:,1);
    TRAJ.gy(j,1).waveform = temp_g_rot(:,2);
end
clear temp_g temp_g_rot;

TRAJ.adc.phaseOffset = pi/2;

%% calculate dummy delay object for reference scans
TRAJ.ref_delay = mr.makeDelay(mr.calcDuration(TRAJ.gx(1), TRAJ.gy(1)));

%% calculate objects for slice excitation
[TRAJ.rf, temp_exc] = mr.makeSincPulse( TRAJ.exc_fa, ...
                                                 system, ...
                                                 'Duration', TRAJ.exc_time,...
                                                 'SliceThickness', TRAJ.slice_thickness, ...
                                                 'maxSlew', system.maxSlew*TRAJ.lim_slew, ...
                                                 'timeBwProduct', TRAJ.exc_tbw, ...
                                                 'apodization', 0.5, ...
                                                 'PhaseOffset', 0, ...
										         'use', 'excitation' );
temp_reph = mr.makeTrapezoid('z', 'Area', -temp_exc.area/2, 'maxGrad', system.maxGrad*TRAJ.lim_grad, 'maxSlew', system.maxSlew*TRAJ.lim_slew, 'system', system);
TRAJ.rf.freqOffset =  TRAJ.slice_offset * temp_exc.amplitude;

%% set channels for gradients
TRAJ.g_exc_x  = temp_exc;  TRAJ.g_exc_x.channel  = 'x';
TRAJ.g_exc_y  = temp_exc;  TRAJ.g_exc_y.channel  = 'y';
TRAJ.g_reph_x = temp_reph; TRAJ.g_reph_x.channel = 'x';
TRAJ.g_reph_y = temp_reph; TRAJ.g_reph_y.channel = 'y';
clear temp_exc temp_reph;

%% calculate spoiler gradient
TRAJ.g_spoil_z = mr.makeTrapezoid('z', 'Area', TRAJ.spoil_nTwists / 1e-3, 'maxGrad', system.maxGrad*TRAJ.lim_grad, 'maxSlew', system.maxSlew*TRAJ.lim_slew, 'system', system);

%% recovery time between repetitions
TRAJ.Trec = mr.makeDelay(round(TRAJ.Trec/system.gradRasterTime)*system.gradRasterTime);

%% delay between objects
TRAJ.d1 = mr.makeDelay(system.gradRasterTime);

%% create dummy FOV object
FOV.Nxy      = 64;
FOV.fov_xy   = 2*TRAJ.slice_offset;
FOV.dz       = TRAJ.slice_thickness;
FOV.z_offset = 0;
FOV_init();

end

