function [a_ramp] = SPI_calc_rf_ramp(NR, a_start)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

    if nargin==2
        [~, a_ramp] = calc_ramp(NR, a_start);
        if max(a_ramp)>pi/2
            warning('alpha ramp error!!! >pi/2')
        end
        if min(a_ramp)==0
            warning('alpha ramp error!!! 0')
        end    
    end
    
    %% last alpha is pi/2
    if (nargin==1 || isempty(a_start))
        a_diff      = @(P) abs( pi/2 - calc_ramp(NR,P(1)) );            
        a_first     = fminsearch(a_diff, 0.1);
        [~, a_ramp] = calc_ramp(NR, a_first);
        Mz_temp     = 1;
        Mz_ramp     = zeros(NR,1);
        Mxy_ramp    = zeros(NR,1);
        for k=1:NR
            Mz_ramp(k,1)  = Mz_temp;
            Mxy_ramp(k,1) = Mz_temp * sin(a_ramp(k));
            Mz_temp = Mz_temp * cos(a_ramp(k));
        end
    end
    
    %%
    function [a_last, a_ramp] = calc_ramp(NR, a_start)
        a_ramp   = zeros(NR,1);
        Mz  = 1;
        Mxy = 0;
        for j=1:NR  
            if j==1
                a = a_start;
            else
                a = asin(Mxy/Mz);
            end
            Mxy         = Mz * sin(a);
            Mz          = Mz * cos(a);
            a_ramp(j)   = a;
        end
        if abs(imag(a_ramp(end)))>0
            clear a_ramp;
            a_ramp = zeros(NR,1);
        end
        a_last = a_ramp(NR);
    end

end

