function [rf] = AHP_pulse( system, mod, tau, wmax, ratio, beta1, beta2, phase_offset )

    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
    
    %  shape based on: adiabatic half passages (AHP)
    %  w(t)  = a * Sech[beta*(t-tau)] + c
    %  dw(t) = a * Tanh[beta*(t-tau/2)] + c
    %  mod: 'tipdown' or 'tipup'
    
    % calculate pulse shape
    n  = round(tau/system.rfRasterTime);
    t  = (1:n)' * system.rfRasterTime;
    f1 = AHP_modulation(t, tau, wmax, ratio, beta1, beta2, mod)/2/pi;
    
    % create dummy pulse object
    rf = mr.makeSincPulse( pi/2, ...
                           system, ...
                           'Duration', tau,...
                           'SliceThickness', 1e6, ...
                           'PhaseOffset', 0, ...
                           'use', 'excitation' );
    
    % replace pulse shape and phase offset
    rf.signal      = rf.signal*0;
    rf.signal(1:n) = f1(:);
    rf.phaseOffset = phase_offset;

end