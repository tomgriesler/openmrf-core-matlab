function [P_corr] = mg_unwrap_phase(S_meas, TEs, unwrap_flag)

% Version: Maximilian Gram, 21.03.2024

    % input
    % S_meas:      complex MRI signal
    % TEs:         echo times
    % unwrap_flag: [0 0] -> no unwrapping
    %              [1 0] -> matlab unwrapping
    %              [0 1] -> mgram unwrapping
    %              [1 1] -> matlab & mgram unwrapping
    % output:
    % P_corr: corrected unwrapped phases
    
    P_corr = angle(S_meas);     
    P_ref  = P_corr(1);
    P_corr = P_corr - P_ref;
    for j=1:numel(P_corr)
        if P_corr(j) > pi
            P_corr(j) = P_corr(j) - 2*pi;
        end
        if P_corr(j) < -pi
            P_corr(j) = P_corr(j) + 2*pi;
        end        
    end
    
    % built-in matlab unwrapper for linear TEs
    if unwrap_flag(1) == 1
        P_corr = double(unwrap(P_corr));
    else
        P_corr = double(P_corr);
    end

    % mgram phase unwrapper for arbitrary TEs
    if unwrap_flag(2) == 1
        TEs    = TEs - TEs(1);
        dw_1st = P_corr(2) / TEs(2);          % 1st order approximation for dw
        P_1st  = dw_1st * TEs;                % 1st order prediction for phase evolution
        n_max  = ceil(abs(dw_1st) * TEs(end) /2/pi);  % maximum number of possible phase wraps
        for n = 1:n_max+1
            for j = 3:numel(TEs)
                if abs(P_1st(j)-P_corr(j)) > pi
                    P_corr(j) = P_corr(j) + 2*pi * sign(dw_1st);  % +- 2pi shift
                end
            end
        end
    end

    % re-add reference phase
    P_corr = P_corr + P_ref;

end

