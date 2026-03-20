function [EPI, ktraj_adc, ktraj_full, ktraj_reco] = EPI_init(EPI, FOV, system)

% ---------------------------------------------------------
% ---------- init parameters and pulseq objects -----------
% -------------- readout: echo-planar (EPI) ---------------
% ---------------------------------------------------------

% source: https://pulseq.github.io/writeEpiRS.html
% this is an experimentaal high-performance EPI sequence
% which uses split gradients to overlap blips with the readout
% gradients combined with ramp-samping
% version 20231112
% M. Gram: deleted slice loop
%        option for SLR pulses

pe_enable         = EPI.pe_enable;
ro_os             = EPI.ro_os;
readoutTime       = EPI.readoutTime;
partFourierFactor = EPI.partFourierFactor;

%% Create 90 degree slice selection pulse and gradient

if strcmp(EPI.exc_mode, 'sinc')
[rf, gz, gzReph] = mr.makeSincPulse( pi/2, 'system', system, ...
                                     'Duration', EPI.exc_time, ...
                                     'SliceThickness', FOV.dz, ...
                                     'apodization', 0.5,...
                                     'timeBwProduct', EPI.exc_tbw, ...,
                                     'apodization', 0.5, ...
                                     'phaseOffset', -pi/2, ...
									 'use', 'excitation' );
 elseif strcmp(EPI.exc_mode, 'sigpy_SLR')
        [rf, gz, gzReph] = SIGPY_SLR( pi/2, EPI.exc_time, -pi/2, EPI.exc_tbw, 'ex', 'ls', 0.01, 0.01, 0 , FOV.dz, system);
end

rf.freqOffset = gz.amplitude * FOV.z_offset;

%% Define other gradients and ADC events
deltak = 1 / FOV.fov_xy;
kWidth = FOV.Nxy * deltak;

%% Phase blip in shortest possible time
blip_dur = ceil(2*sqrt(deltak/system.maxSlew)/system.gradRasterTime/2) * system.gradRasterTime * 2;    % we round-up the duration to 2x the gradient raster time
gy       = mr.makeTrapezoid('y', system, 'Area', -deltak, 'Duration', blip_dur);                       % we use negative blips to save one k-space line on our way towards the k-space center

%% readout gradient is a truncated trapezoid with dead times at the beginnig
% and at the end each equal to a half of blip_dur
% the area between the blips should be defined by kWidth
% we do a two-step calculation: we first increase the area assuming maximum
% slewrate and then scale down the amlitude to fix the area
extra_area = blip_dur / 2 * blip_dur / 2 *system.maxSlew;  % check unit!;
gx           = mr.makeTrapezoid('x', system, 'Area', kWidth+extra_area, 'duration', readoutTime+blip_dur);
actual_area  = gx.area-gx.amplitude / gx.riseTime * blip_dur / 2 * blip_dur/2/2 - gx.amplitude / gx.fallTime * blip_dur / 2 * blip_dur/2/2;
gx.amplitude = gx.amplitude / actual_area*kWidth;
gx.area      = gx.amplitude * (gx.flatTime + gx.riseTime/2 + gx.fallTime/2);
gx.flatArea  = gx.amplitude * gx.flatTime;

%% calculate ADC
% we use ramp sampling, so we have to calculate the dwell time and the
% number of samples, which are will be qite different from Nx and
% readoutTime/Nx, respectively.
adcDwellNyquist = deltak / gx.amplitude / ro_os;
adcDwell        = floor(adcDwellNyquist*1e7) * 1e-7;
adcSamples      = floor(readoutTime/adcDwell/4) *4; % on Siemens the number of ADC samples need to be divisible by 4
adc             = mr.makeAdc(adcSamples, 'Dwell', adcDwell, 'Delay', blip_dur/2);
time_to_center  = adc.dwell*((adcSamples-1)/2+0.5); % I've been told that Siemens samples in the center of the dwell period
adc.delay       = round((gx.riseTime+gx.flatTime/2-time_to_center)*1e6)*1e-6; % we adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us
% this rounding actually makes the sampling points on odd and even readouts
% to appear misalligned. However, on the real hardware this misalignment is
% much stronger anyways due to the grdient delays

%% split the blip into two halves and produnce a combined synthetic gradient
gy_parts                 = mr.splitGradientAt(gy, blip_dur/2, system);
[gy_blipup, gy_blipdown] = mr.align('right',gy_parts(1), 'left', gy_parts(2), gx);
gy_blipdownup            = mr.addGradients({gy_blipdown, gy_blipup}, system);

%% pe_enable support
gy_blipup.waveform     = gy_blipup.waveform     * pe_enable;
gy_blipdown.waveform   = gy_blipdown.waveform   * pe_enable;
gy_blipdownup.waveform = gy_blipdownup.waveform * pe_enable;

%% phase encoding and partial Fourier
Ny_pre  = round(partFourierFactor*FOV.Nxy/2-1); % PE steps prior to ky = 0, excluding the central line
Ny_post = round(FOV.Nxy/2+1); % PE lines after the k-space center including the central line
Ny_meas = Ny_pre + Ny_post;

%% Pre-phasing gradients
gxPre                  = mr.makeTrapezoid('x', system, 'Area', -gx.area/2);
gyPre                  = mr.makeTrapezoid('y', system, 'Area', Ny_pre*deltak);
[gxPre, gyPre, gzReph] = mr.align('right', gxPre, 'left', gyPre,gzReph);
gyPre                  = mr.makeTrapezoid('y', system, 'Area', gyPre.area, 'Duration', mr.calcDuration(gxPre,gyPre,gzReph));
gyPre.amplitude        = gyPre.amplitude * pe_enable;

%% save EPI ojects
EPI.rf            = rf;
EPI.gz            = gz;
EPI.gxPre         = gxPre;
EPI.gyPre         = gyPre;
EPI.gzReph        = gzReph;
EPI.gx            = gx;
EPI.gy_blipup     = gy_blipup;
EPI.gy_blipdown   = gy_blipdown;
EPI.gy_blipdownup = gy_blipdownup;
EPI.adc           = adc;
EPI.Ny_meas       = Ny_meas;

%% calculate kspace trajectory
seq = mr.Sequence(system);
EPI_add();
[ktraj_adc, ktraj_full] = pulseq_get_ktraj(seq, 1);
ktraj_x                 = ktraj_adc(1,:)';
ktraj_y                 = ktraj_adc(2,:)';
ktraj_x                 = reshape(ktraj_x, EPI.adc.numSamples, EPI.Ny_meas)';
ktraj_y                 = reshape(ktraj_y, EPI.adc.numSamples, EPI.Ny_meas)';
ktraj_reco(1,:,:)       = ktraj_x(:,:);
ktraj_reco(2,:,:)       = ktraj_y(:,:);

end
