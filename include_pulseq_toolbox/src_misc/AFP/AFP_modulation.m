function w1 = AFP_modulation(t, tau, wmax, ratio, beta1, beta2)

    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

    % re-scale inputs   
    if size(t,1)<size(t,2)
        t = t';
    end
    dwmax = wmax * ratio;
    b1    = beta1 * 1e3;
    b2    = beta2 * 1e3;

    % modulations for w1 profile: sech()
    w1_1 = wmax+(wmax-wmax*sech(b1*(t-tau)))/(-1+sech(b1*tau));  % tipdown  
    w1_2 = wmax+(wmax-wmax*sech(b1*t))/(-1+sech(b1*tau));        % tipup
    w1   = [w1_1; w1_2];
    
    % modulations for phase profile: tanh()
    phi_1 = (dwmax*(b2*(t-tau)+coth((b2*tau)/2)*(-log(cosh(b2*(t-tau/2)))+log(cosh((b2*tau)/2)))))/(2*b2);   % tipdown
    phi_2 = (dwmax*(b2*t+coth((b2*tau)/2)*(log(cosh(b2*(t-tau/2)))-log(cosh((b2*tau)/2)))))/(2*b2);          % tipup
    phi_1 = phi_1 - phi_1(end);
    phi_2 = phi_2 - phi_2(1);
    phi_2 = -phi_2;
    phi   = [phi_1; phi_2];

    % output complex w1 waveform
    w1 = w1 .* exp(1i*phi);

end

