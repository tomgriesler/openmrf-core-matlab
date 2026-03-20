function [AHP] = AHP_get_params(tau, wmax)
% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
    ahp_path = pulseq_get_path('AHP_get_params');
    ahp_path = [ahp_path 'AHP_params_' num2str(tau*1e3) 'ms_' num2str(wmax/2/pi) 'Hz.mat'];
    if exist(ahp_path, 'file')
        load(ahp_path, 'ratio', 'beta1', 'beta2');
    else
        P     = AHP_optimize_params(tau, wmax);
        ratio = P(1);
        beta1 = P(2);
        beta2 = P(3);
        save(ahp_path, 'ratio', 'beta1', 'beta2');
    end
    AHP.tau   = tau;
    AHP.wmax  = wmax;
    AHP.ratio = ratio;
    AHP.dwmax = wmax * ratio;
    AHP.beta1 = beta1;
    AHP.beta2 = beta2;
end

