function [rf, gz, gz_reph] = SIGPY_SLR(flip_angle, pulse_duration, phase_offset, tbw, ptype, ftype, d1, d2, cancel_alpha_phs, dz, sys, maxSlew)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% calculate pulse object using: sigpy.mri.rf.slr.dzrf
% sigpy-parameters:
%    n       number of time points
%    tbw     pulse time bandwidth product
%    ptype   pulse type:  st  (small-tip excitation)
%                         ex  (pi/2 excitation pulse)
%                         se  (spin-echo pulse)
%                         inv  (inversion)
%                         sat  (pi/2 saturation pulse)
%    ftype   type of filter to use:  ms  (sinc)
%                                    pm  (Parks-McClellan equal-ripple)
%                                    min  (minphase using factored pm)
%                                    max  (maxphase using factored pm)
%                                    ls  (least squares)
%    d1      passband ripple level
%    d2      stopband ripple level
%    cancel_alpha_phs - bool
%    for  ex  pulses, absorb the alpha phase profile from beta s profile,
%    so they cancel for a flatter total phase

if nargin<12
    maxSlew = [];
end
if isempty(maxSlew)
    maxSlew = sys.maxSlew;
end
sys.maxSlew = maxSlew;

n = round(pulse_duration/sys.rfRasterTime);

%% call python-clone function: only works for ftype = 'ls'
signal = SIGPY_slr_dzrf(n, tbw, ptype, ftype, d1, d2, cancel_alpha_phs);

%% create sigpy pulse object and gradients

% scale waveform for flip angle
signal = signal / (abs(sum(signal)) * sys.rfRasterTime *2*pi) *flip_angle;

% create dummy pulse object
[rf, gz, gz_reph] = mr.makeSincPulse( pi/100, ... % only a dummy
                                      sys, ...
                                      'Duration',       pulse_duration,...
                                      'SliceThickness', dz, ...
                                      'timeBwProduct',  tbw, ...
                                      'PhaseOffset',    0, ...
                                      'use', 'excitation' );

% replace pulse shape and phase offset
rf.signal      = rf.signal*0;
rf.signal(1:n) = signal(:);
rf.phaseOffset = phase_offset;

% rf use
if strcmp(ptype, 'st')
    rf.use = 'excitation';
end
if strcmp(ptype, 'ex')
    rf.use = 'excitation';
end
if strcmp(ptype, 'se')
    rf.use = 'refocusing';
end
if strcmp(ptype, 'inv')
    rf.use = 'inversion';
end
if strcmp(ptype, 'sat')
    rf.use = 'saturation';
end    

end