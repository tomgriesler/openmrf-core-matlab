function [rf] = SIGPY_HYPSEC(beta, mu, pulse_duration, phase_offset, sys)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% calculate pulse object using: sigpy.mri.rf.adiabatic.hypsec
% sigpy-parameters:
%    n       number of time points
%    beta  - AM waveform parameter
%    mu    - a constant, determines amplitude of frequency sweep

n = round(pulse_duration/sys.rfRasterTime);

%% calculate sigpy pulse shape

[signal_am, signal_fm] = hypsec_clone(n, beta, mu, pulse_duration);
signal_am = signal_am / max(signal_am) * beta; 

% calculate phase modulation
signal_phase = cumsum(signal_fm) * sys.rfRasterTime;

% calculate complex pulse shape
signal = signal_am .* exp(1i* signal_phase);
       


%% create sigpy pulse object

% create dummy pulse object
rf = mr.makeSincPulse( pi, ...
                       sys, ...
                       'Duration', pulse_duration,...
                       'SliceThickness', 1e6, ...
                       'PhaseOffset', 0, ...
                       'use', 'inversion' );

% replace pulse shape and phase offset
rf.signal      = rf.signal*0;
rf.signal(1:n) = signal(:);
rf.phaseOffset = phase_offset;

end

%% clone of the python based SIGPY hypsech function
function [a, om] = hypsec_clone(n, beta, mu, dur)
% Design a hyperbolic secant adiabatic pulse.
%
% mu * beta becomes the amplitude of the frequency sweep
%
% Args:
%     n (int): number of samples (should be a multiple of 4).
%     beta (float): AM waveform parameter.
%     mu (float): a constant, determines amplitude of frequency sweep.
%     dur (float): pulse time (s).
%
% Returns:
%     2-element tuple containing
%      - a (array): AM waveform.
%      - om (array): FM waveform (radians/s).
%
% References:
%     Baum, J., Tycko, R. and Pines, A. (1985). 'Broadband and adiabatic
%     inversion of a two-level system by phase-modulated pulses'.
%     Phys. Rev. A., 32:3435-3447.

t = ((-n/2):(n/2-1)).' / n * dur;
a = (cosh(beta * t)).^(-1);
om = -mu * beta * tanh(beta * t);
end
