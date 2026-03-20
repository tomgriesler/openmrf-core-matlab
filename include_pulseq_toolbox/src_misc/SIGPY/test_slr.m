%% compare: sigpy SLR vs pulseq sinc
clear
pulseq_init();
flip_angle     = 20 *pi/180;
pulse_duration = 1 *1e-3;
tb             = 6;
ptype          = 'st';
dz             = 8 *1e-3;

phase_offset   = 0;
ftype          = 'ls';
d1             = 0.01;
d2             = 0.01;
cap            = 0;

[rf_sigpy,  gz_sigpy]  = SIGPY_SLR(flip_angle, pulse_duration, phase_offset, tb, ptype, ftype, d1, d2, cap ,dz, system);
[rf_pulseq, gz_pulseq] = mr.makeSincPulse( flip_angle, ...
                         system, ...
                         'Duration', pulse_duration, ...
                         'SliceThickness', dz, ...
                         'timeBwProduct', tb, ...
                         'PhaseOffset', phase_offset, ...
                         'apodization', 0.5, ...
						 'use', 'excitation' );

%%
figure()
hold on
plot(rf_pulseq.t*1e3, abs(rf_pulseq.signal), 'LineWidth', 2)
plot(rf_sigpy.t*1e3,  abs(rf_sigpy.signal),  'LineWidth', 2)
legend('pulseq sinc', 'sigpy slr')
xlabel('time t [ms]')
ylabel('abs(f1) [Hz]')
hold off

%% show slice simulation
mg_slice_sim(rf_pulseq.signal, dz, gz_pulseq.amplitude, flip_angle, system.rfRasterTime, 1000, 2);
mg_slice_sim(rf_sigpy.signal,  dz, gz_sigpy.amplitude,  flip_angle, system.rfRasterTime, 1000, 2);
