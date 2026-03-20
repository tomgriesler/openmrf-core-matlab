function FAs = MRF_calc_FAs_sin(FA_list)

    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

    % ----- input: -----
    % FA_list: Nx3 array, [min1, max1, n1; min2, max2, n2; ...]
    
    % ----- output: -----
    % FAs: N*n x 1 array, sinusoidal flip anlge pattern
    
    FAs = [];
    for j = 1:size(FA_list,1)
        FAs = [FAs; sin(linspace(0,pi,FA_list(j,3))') * (FA_list(j,2)-FA_list(j,1)) + FA_list(j,1)];
    end

end