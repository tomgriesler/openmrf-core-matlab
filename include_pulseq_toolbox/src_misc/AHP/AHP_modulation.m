function [w1, u1] = AHP_modulation(t, tau, wmax, ratio, beta1, beta2, mod)

    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

    % re-scale inputs   
    if size(t,1)<size(t,2)
        t = t';
    end
    dwmax = wmax * ratio;
    b1    = beta1 * 1e3;
    b2    = beta2 * 1e3;
    
    % modulations for w1 profile: sech()    
    if strcmp(mod, 'tipdown')
        w1 = wmax+(wmax-wmax*sech(b1*(t-tau)))/(-1+sech(b1*tau));  
    elseif strcmp(mod, 'tipup')
        w1 = wmax+(wmax-wmax*sech(b1*t))/(-1+sech(b1*tau));  
    end
    
    % modulations for phase profile: tanh()
    if strcmp(mod, 'tipdown')
        phi = (dwmax*(b2*(t-tau)+coth((b2*tau)/2)*(-log(cosh(b2*(t-tau/2)))+log(cosh((b2*tau)/2)))))/(2*b2);
    elseif strcmp(mod, 'tipup')
        phi = (dwmax*(b2*t+coth((b2*tau)/2)*(log(cosh(b2*(t-tau/2)))-log(cosh((b2*tau)/2)))))/(2*b2);
    end
    
    % modulations for offresonance profile: tanh()
    if strcmp(mod, 'tipdown')
        dw = 1/2*(dwmax-dwmax*coth((b2*tau)/2)*tanh(b2*(t-tau/2)));
    elseif strcmp(mod, 'tipup')
        dw = 1/2*(dwmax+dwmax*coth((b2*tau)/2)*tanh(b2*(t-tau/2)));
    end
       
    % calculate unit vector (direction) of effective field
    u1(1,:) = w1 .* cos(phi);
    u1(2,:) = w1 .* sin(phi);
    u1(3,:) = dw;
    u1      = u1 ./ sqrt( u1(1,:).^2 + u1(2,:).^2 + u1(3,:).^2 );

    % output complex w1 waveform
    w1 = w1 .* exp(1i*phi);

end