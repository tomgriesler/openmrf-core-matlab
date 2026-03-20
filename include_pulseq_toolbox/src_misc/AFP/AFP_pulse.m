function [rf] = AFP_pulse( system, tau, wmax, ratio, beta1, beta2, phase_offset )

    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
    
    %  shape based on: adiabatic half passages (AHP)
    %  w(t)  = a * Sech[beta*(t-tau)] + c
    %  dw(t) = a * Tanh[beta*(t-tau/2)] + c
    %  use ahp tipdown for 1st half
    %  use ahp tipup for 2nd half, but with inversed phase modulation
    
    % calculate pulse shape
    n  = round(tau/system.rfRasterTime);
    t  = (1:n)' * system.rfRasterTime;
    f1 = AFP_modulation(t, tau, wmax, ratio, beta1, beta2)/2/pi;
    n  = numel(f1);
    
    % create dummy pulse object
    rf = mr.makeSincPulse( pi, ...
                           system, ...
                           'Duration', 2*tau,...
                           'SliceThickness', 1e6, ...
                           'PhaseOffset', 0, ...
                           'use', 'inversion' );
    
    % replace pulse shape and phase offset
    rf.signal      = rf.signal*0;
    rf.signal(1:n) = f1(:);
    rf.phaseOffset = phase_offset;

end