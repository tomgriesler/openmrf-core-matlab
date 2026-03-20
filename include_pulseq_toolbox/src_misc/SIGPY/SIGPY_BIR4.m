function [rf] = SIGPY_BIR4(beta, kappa, dw0, f1, alpha, dphi, pulse_duration, phase_offset, mod, sys)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% calculate pulse object using: sigpy.mri.rf.adiabatic.bir4
%    n       number of time points
%    beta    AM waveform parameter
%    kappa - FM waveform parameter
%    alpha   flip angle [rad]
%    dphi  - phase tuning [rad]
%    dw0     maximum offresonance [rad/s]
%    f1    - maximum B1+ field strenth [Hz]
%    mod   - 'tipdown' or 'tipup'

n = round(pulse_duration/sys.rfRasterTime);

%% calculate sigpy pulse shape

[signal_am, signal_fm] = bir4_clone(n, beta, kappa, alpha, dw0);
signal_am = signal_am / max(signal_am) * f1;

% calculate phase modulation
signal_phase = cumsum(signal_fm) * sys.rfRasterTime;

% calculate complex pulse shape
if strcmp(mod, 'tipdown')
    signal = signal_am .* exp(1i* (signal_phase + dphi - pi/2));
elseif strcmp(mod, 'tipup')
    signal = signal_am .* exp(1i* (signal_phase + dphi + pi/2));
else
    error('wrong mode: tipdown or tipup')
end

%% create sigpy pulse object

% create dummy pulse object
rf = mr.makeSincPulse( pi/2, ...
                       sys, ...
                       'Duration', pulse_duration,...
                       'SliceThickness', 1e6, ...
                       'PhaseOffset', 0, ...
                       'use', 'excitation' );

% replace pulse shape and phase offset
rf.signal      = rf.signal*0;
rf.signal(1:n) = signal(:);
rf.phaseOffset = phase_offset;

end

%% clone of the python based SIGPY bir4 function
function [a, om] = bir4_clone(n, beta, kappa, theta, dw0)
% Design a BIR-4 adiabatic pulse.
%
% BIR-4 is equivalent to two BIR-1 pulses back-to-back.
%
% Args:
%     n (int): number of samples (should be a multiple of 4).
%     beta (float): AM waveform parameter.
%     kappa (float): FM waveform parameter.
%     theta (float): flip angle in radians.
%     dw0: FM waveform scaling (radians/s).
%
% Returns:
%     2-element tuple containing
%      - a (array): AM waveform.
%      - om (array): FM waveform (radians/s).
%
% References:
%     Staewen, R.S. et al. (1990). '3-D FLASH Imaging using a single surface
%     coil and a new adiabatic pulse, BIR-4'.
%     Invest. Radiology, 25:559-567.

dphi = pi + theta / 2;
t = (0:n-1).' / n;

q = n/4;
a1 = tanh(beta * (1 - 4 * t(1:q)));
a2 = tanh(beta * (4 * t(q+1:2*q) - 1));
a3 = tanh(beta * (3 - 4 * t(2*q+1:3*q)));
a4 = tanh(beta * (4 * t(3*q+1:n) - 3));

a = [a1; a2; a3; a4];
a = complex(a, 0);
a(q+1:3*q) = a(q+1:3*q) .* exp(1i * dphi);

om1 = dw0 * tan(kappa * 4 * t(1:q)) / tan(kappa);
om2 = dw0 * tan(kappa * (4 * t(q+1:2*q) - 2)) / tan(kappa);
om3 = dw0 * tan(kappa * (4 * t(2*q+1:3*q) - 2)) / tan(kappa);
om4 = dw0 * tan(kappa * (4 * t(3*q+1:n) - 4)) / tan(kappa);

om = [om1; om2; om3; om4];
end
