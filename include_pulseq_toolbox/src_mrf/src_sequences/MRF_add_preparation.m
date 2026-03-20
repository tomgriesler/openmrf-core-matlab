% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% MRF preparations

% ----- No Preparation -----
if strcmp(MRF.enc_list{loop_MRF}, 'No_Prep')
    % no preparation
end

% ----- Saturation -----
if strcmp(MRF.enc_list{loop_MRF}, 'Saturation')
    SAT_add();
    loop_SAT = loop_SAT + 1;
end

% ----- Inversion -----
if strcmp(MRF.enc_list{loop_MRF}, 'Inversion')
    INV_add();
    loop_INV = loop_INV + 1;
end

% ----- T2 Preparation -----
if strcmp(MRF.enc_list{loop_MRF}, 'T2')
    T2_add();
    loop_T2 = loop_T2 + 1;
end

% ----- Spin-Lock Preparation -----
if strcmp(MRF.enc_list{loop_MRF}, 'SL')
    SL_add();
    loop_SL = loop_SL + 1;
end

% ----- MLEV Preparation -----
if strcmp(MRF.enc_list{loop_MRF}, 'MLEV')
    MLEV_add();
    loop_MLEV = loop_MLEV + 1;
end

% ----- adiabatic Spin-Lock Preparation -----
if strcmp(MRF.enc_list{loop_MRF}, 'ADIASL')
    ADIASL_add();
    loop_ADIASL = loop_ADIASL + 1;
end    

% ----- Fat Saturation -----
if exist('FAT', 'var')
    if ( strcmp(MRF.enc_list{loop_MRF}, 'Inversion') && INV.inv_rec_time(loop_INV-1) < 0.2 )
        % no fat saturation
    else    
        FAT_add();
    end
end