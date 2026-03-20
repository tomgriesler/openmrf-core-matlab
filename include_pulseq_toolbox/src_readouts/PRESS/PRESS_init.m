function PRESS = PRESS_init(PRESS, FOV, system, flag_plot)

% ---------------------------------------------------------
% ---------- init parameters and pulseq objects -----------
% ----- readout: Point RESolved Spectroscopy (PRESS) ------
% ---------------------------------------------------------

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

if nargin < 4
    flag_plot = 0;
end

%% PRESS calculate RF objects

if strcmp(PRESS.exc.mode, 'sinc')

[PRESS.exc.rf90z, PRESS.exc.g90z] = mr.makeSincPulse( pi/2, ...
                                                      system, ...
                                                      'Duration', PRESS.exc.t90, ...
                                                      'timeBwProduct', PRESS.exc.tbw, ...
                                                      'PhaseOffset', PRESS.exc.phz, ...
                                                      'SliceThickness', FOV.dz, ...
                                                      'maxSlew', system.maxSlew*0.9, ...
                                                      'use', 'excitation');
PRESS.exc.rf90z.freqOffset = PRESS.exc.g90z.amplitude * FOV.z_offset;
PRESS.exc.g90z.channel     = 'z';
PRESS.exc.g90z_reph        = mr.makeTrapezoid('z', 'Area', -PRESS.exc.g90z.area/2, 'Duration', PRESS.exc.t90/2, 'maxGrad', system.maxGrad*0.8, 'maxSlew', system.maxSlew*0.8, 'system', system);

[PRESS.exc.rf180y, PRESS.exc.g180y] = mr.makeSincPulse( pi, ...
                                                        system, ...
                                                        'Duration', PRESS.exc.t180, ...
                                                        'timeBwProduct', PRESS.exc.tbw, ...
                                                        'PhaseOffset', PRESS.exc.phy, ...
                                                        'SliceThickness', FOV.dy, ...
                                                        'maxSlew', system.maxSlew*0.9, ...
                                                        'use', 'refocusing');
PRESS.exc.rf180y.freqOffset = PRESS.exc.g180y.amplitude * FOV.y_offset;
PRESS.exc.g180y.channel    = 'y';

[PRESS.exc.rf180x, PRESS.exc.g180x] = mr.makeSincPulse( pi, ...
                                                        system, ...
                                                        'Duration', PRESS.exc.t180, ...
                                                        'timeBwProduct', PRESS.exc.tbw, ...
                                                        'PhaseOffset', PRESS.exc.phx, ...
                                                        'SliceThickness', FOV.dx, ...
                                                        'maxSlew', system.maxSlew*0.9, ...
                                                        'use', 'refocusing');
PRESS.exc.rf180x.freqOffset = PRESS.exc.g180x.amplitude * FOV.x_offset;
PRESS.exc.g180x.channel     = 'x';
    
elseif strcmp(GRE.exc_mode, 'sigpy_SLR')

error('work in progress!')

end

%% PRESS calculate CRUSH objects
PRESS.crush.gy = mr.makeTrapezoid('y', 'Area', -PRESS.crush.n_twists / FOV.dy, 'Duration', PRESS.crush.duration, 'maxGrad', system.maxGrad*0.8, 'maxSlew', system.maxSlew*0.8, 'system', system);
PRESS.crush.gx = mr.makeTrapezoid('x', 'Area', -PRESS.crush.n_twists / FOV.dx, 'Duration', PRESS.crush.duration, 'maxGrad', system.maxGrad*0.8, 'maxSlew', system.maxSlew*0.8, 'system', system);

%% PRESS calculate ADC object
if isfield(PRESS.acq, 'adcBWDesired')
    PRESS.acq.adcNSamplesDesired = ceil(PRESS.acq.adcTimeDesired * PRESS.acq.adcBWDesired);
elseif isfield(PRESS.acq, 'adcNSamplesDesired')
    PRESS.acq.adcBWDesired = PRESS.acq.adcNSamplesDesired / PRESS.acq.adcTimeDesired;
end
PRESS.acq.adcTime         = system.gradRasterTime * round( PRESS.acq.adcTimeDesired / system.gradRasterTime ) - system.gradRasterTime;
PRESS.acq.adcNSamples     = ceil(PRESS.acq.adcNSamplesDesired/64)*64;
PRESS.acq.adcDwell        = round( PRESS.acq.adcTime / PRESS.acq.adcNSamples /100e-9 ) *100e-9;
PRESS.acq.adcBW           = 1 / PRESS.acq.adcDwell;
PRESS.acq.adc             = mr.makeAdc(PRESS.acq.adcNSamples, 'Dwell', PRESS.acq.adcDwell); 
PRESS.acq.adc.delay       = system.adcDeadTime;
PRESS.acq.adc.phaseOffset = 0;
PRESS.acq.delay           = round( (mr.calcDuration(PRESS.acq.adc)+1e-3)/system.gradRasterTime ) * system.gradRasterTime;
PRESS.acq.delay           = mr.makeDelay(PRESS.acq.delay);

%% PRESS calculate TE delay
iterate_TE = 1;
while iterate_TE==1
    temp_dt_90_180   = mr.calcDuration(PRESS.exc.rf90z)/2 + mr.calcDuration(PRESS.exc.g90z_reph) + mr.calcDuration(PRESS.crush.gy) + mr.calcDuration(PRESS.exc.rf180y)/2;
    temp_dt_180_180  = mr.calcDuration(PRESS.exc.rf180y)/2 + mr.calcDuration(PRESS.crush.gy) + mr.calcDuration(PRESS.crush.gx) + mr.calcDuration(PRESS.exc.rf180x)/2;
    temp_dt_90_180   = round( (PRESS.TE.echo_filling_delay - temp_dt_90_180) / system.gradRasterTime ) * system.gradRasterTime;
    temp_dt_180_180  = round( (2*PRESS.TE.echo_filling_delay - temp_dt_180_180) / system.gradRasterTime ) * system.gradRasterTime;
    if temp_dt_90_180 > 2*system.gradRasterTime && temp_dt_180_180 > 2*system.gradRasterTime
        iterate_TE = 0;
    else
        PRESS.TE.echo_filling_delay = PRESS.TE.echo_filling_delay + 1e-3;
    end
end
PRESS.TE.delay1 = mr.makeDelay( temp_dt_90_180 );
PRESS.TE.delay2 = mr.makeDelay( temp_dt_180_180 );
clear temp_dt_90_180 temp_dt_180_180;

%% check echo timings

PRESS.TE.dt_90_180 = ...
mr.calcDuration(PRESS.exc.rf90z, PRESS.exc.g90z)/2 + ...
mr.calcDuration(PRESS.exc.g90z_reph) + ...
mr.calcDuration(PRESS.TE.delay1) + ...
mr.calcDuration(PRESS.crush.gy) + ...
mr.calcDuration(PRESS.exc.rf180y, PRESS.exc.g180y)/2;

PRESS.TE.dt_180_180 = ...
mr.calcDuration(PRESS.exc.rf180y, PRESS.exc.g180y)/2 + ...    
mr.calcDuration(PRESS.crush.gy) + ...
mr.calcDuration(PRESS.TE.delay2) + ...
mr.calcDuration(PRESS.crush.gx) + ...
mr.calcDuration(PRESS.exc.rf180x, PRESS.exc.g180x)/2;

if abs(PRESS.TE.dt_180_180 / PRESS.TE.dt_90_180 - 2) > 1e-3
    error('echo timings incorrect!');
end

%% check echo position

if PRESS.TE.echo_pos_delay <= 0
    PRESS.TE.echo_pos_delay = system.gradRasterTime;
end

PRESS.TE.echo_pos_delay   = round(PRESS.TE.echo_pos_delay / system.gradRasterTime) * system.gradRasterTime;
PRESS.TE.TE_eff           = PRESS.TE.dt_90_180 + PRESS.TE.dt_180_180 + PRESS.TE.dt_90_180;
PRESS.TE.delay3           = mr.makeDelay( PRESS.TE.echo_pos_delay );
PRESS.TE.dt_180_echo      = PRESS.TE.dt_90_180;
PRESS.TE.dt_180_adc_start = ...
                            mr.calcDuration(PRESS.exc.rf180x, PRESS.exc.g180x)/2 + ...
                            mr.calcDuration(PRESS.crush.gx) + ...
                            mr.calcDuration(PRESS.TE.delay3);
PRESS.TE.dt_180_adc_end   = ...
                            mr.calcDuration(PRESS.exc.rf180x, PRESS.exc.g180x)/2 + ...
                            mr.calcDuration(PRESS.crush.gx) + ...
                            mr.calcDuration(PRESS.TE.delay3) + ...
                            mr.calcDuration(PRESS.acq.adc);

if PRESS.TE.dt_180_echo < PRESS.TE.dt_180_adc_start
    warning('missing echo in ADC! -> decrease PRESS.TE.echo_pos_delay');
end
if PRESS.TE.dt_180_adc_end < PRESS.TE.dt_180_echo
    warning('missing echo in ADC! -> increase PRESS.TE.echo_pos_delay');
end

PRESS.TE.echo_pos = (PRESS.TE.dt_180_echo - PRESS.TE.dt_180_adc_start) / (PRESS.TE.dt_180_adc_end - PRESS.TE.dt_180_adc_start);

%% show timings
if flag_plot==1
figure()
title(['PRESS timings, echo position: ' num2str(PRESS.TE.echo_pos*100, '%.1f') '%' ])
hold on
xline(0,                                                                      'r-',  'LineWidth', 3, 'Label', '90')
xline(PRESS.TE.dt_90_180 *1e3,                                                'b-',  'LineWidth', 3, 'Label', '180')
xline((PRESS.TE.dt_90_180 + PRESS.TE.dt_180_180) *1e3,                        'b-',  'LineWidth', 3, 'Label', '180')
xline((PRESS.TE.dt_90_180 + PRESS.TE.dt_180_180 + PRESS.TE.dt_180_echo) *1e3, 'r--', 'LineWidth', 3, 'Label', 'Echo')
temp_adc = linspace(PRESS.TE.dt_90_180 + PRESS.TE.dt_180_180 + PRESS.TE.dt_180_adc_start, PRESS.TE.dt_90_180 + PRESS.TE.dt_180_180 + PRESS.TE.dt_180_adc_end, 100) *1e3;
plot(temp_adc, zeros(1,100), 'kx')
xlabel('time [ms]')
ylim([-0.5 0.5])
set(gca,'ytick',[])
clear temp_adc;
end

end