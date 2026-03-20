function [rf, gz] = SIGPY_GOIA_WURST(dur, f, n_b1, m_grad, b1_max, dz, phase_offset, t_rise, t_spoil, system)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% Input:
% dur:          [s]   pulse duration
% f:            [0,1] gradient modulation factor
% n_b1:         [ ]   order for B1 modulation
% m_grad:       [ ]   order for gradient modulation
% b1_max:       [Hz]  maximum B1 amplitude
% dz:           [m]   slice thickness
% phase_offset: [rad] phase offset
% t_rise:       [s]   rise time
% t_spoil:      [s]   duration of each spoiler

bw = system.maxGrad * dz;
n  = round(dur/system.rfRasterTime);

%% modulation functions were copied from SIGPY:
% sigpy.readthedocs.io/en/latest/_modules/sigpy/mri/rf/adiabatic.html#goia_wurst
t  = linspace(0,n,n) * dur/n;
am = b1_max * (1 - abs(sin(pi / 2 * (2 * t / dur - 1))) .^ n_b1);  % [Hz] amplitude modulation
gm = (1 - f) + f * abs(sin(pi / 2 * (2 * t / dur - 1))) .^ m_grad; % [  ] norm. gradient modulation
om = cumsum((am.^2) ./ gm) * dur/n;                                % [Hz] off-resonance modulation
om = om - om(round(n/2)+1);
om = gm .* om;
om = om / max(abs(om)) * bw / 2;
pm = 2*pi * cumsum(om)*system.rfRasterTime; % [rad] phase modulation

% Reference:
% O. C. Andronesi, S. Ramadan, E.-M. Ratai, D. Jennings, C. E. Mountford,
% A. G. Sorenson.
% J Magn Reson, 203:283-293, 2010.

%% create rf and gz pulseq objects

% create dummy pulse object
rf = mr.makeSincPulse( pi, ...
                       system, ...
                       'Duration', dur,...
                       'PhaseOffset', phase_offset, ...
                       'use', 'inversion' );

% replace amplitude and phase modulation
rf.signal = am .* exp(1i*pm);

% scaling gradient waveform
gmax = system.maxGrad;
gz   = [linspace(0, 1, round(t_rise/1e-5)), ...
          ones(1,round(t_spoil/1e-5)), ...
          interp1(linspace(0,1,n), gm, linspace(0,1,round(n*system.rfRasterTime/system.gradRasterTime))), ...
          ones(1,round(t_spoil/1e-5)), ...
          linspace(1, 0, round(t_rise/1e-5))] * gmax;
gz = mr.makeArbitraryGrad('z', gz, 'system', system, 'first', 0, 'last', 0);
rf.delay = t_rise + t_spoil;

end