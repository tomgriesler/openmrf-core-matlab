function FAs = MRF_calc_FAs_sin_rand(FA_min, FA_max, nr, n_segm, f)

    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

    if nargin<5
        f = 0.7;
    end

    rng("default");    
    FAs = zeros(nr, n_segm);
    for j = 1:n_segm
        temp_max = (FA_max-FA_min)*2/3 * rand() + (FA_max-FA_min)/3 + FA_min;
        temp_fa  = FA_min + (temp_max-FA_min) * sin(linspace(0, pi*f, nr));
        FAs(:,j) = temp_fa(:);
    end
    FAs = FAs(:);


end

