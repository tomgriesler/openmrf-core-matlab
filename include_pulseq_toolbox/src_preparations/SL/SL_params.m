% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% standard params for SL preparation

SL.freqOffset = 0.0; % [Hz] frequency offset of SL pulses

% excitation pulses
SL.exc_mode  = 'adiabatic_AHP';  % 'adiabatic_AHP', 'adiabatic_BIR4', 'sinc', 'sigpy_SLR' or 'bp'
SL.exc_time  = 3.0 *1e-3;        % [s] excitation time
SL.exc_tbw   = [];               % [ ] time bandwith product, only for sinc/sigpy_SLR
SL.adia_wmax = 600 * 2*pi;       % [rad/s] amplitude of adiabatic pulse

% refocusing pulses
SL.rfc_mode  = 'sinc';           % 'sinc', 'sigpy_SLR', 'bp' or 'comp'
SL.rfc_time  = 4.0 *1e-3;        % [s] refocusing time
SL.rfc_tbw   = 4;                % [ ] time bandwith product, only for sinc/sigpy_SLR
