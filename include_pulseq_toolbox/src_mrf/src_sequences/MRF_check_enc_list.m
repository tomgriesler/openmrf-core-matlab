% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% check MRF encoding list

MRF.n_no_prep   = sum(strcmp(MRF.enc_list, 'No_Prep'));
MRF.n_sat       = sum(strcmp(MRF.enc_list, 'Saturation'));
MRF.n_inversion = sum(strcmp(MRF.enc_list, 'Inversion'));
MRF.n_sl        = sum(strcmp(MRF.enc_list, 'SL'));
MRF.n_mlev      = sum(strcmp(MRF.enc_list, 'MLEV'));
MRF.n_t2        = sum(strcmp(MRF.enc_list, 'T2'));
MRF.n_adiasl    = sum(strcmp(MRF.enc_list, 'ADIASL'));

if MRF.n_segm ~= MRF.n_no_prep + MRF.n_sat + MRF.n_inversion + MRF.n_sl + MRF.n_mlev + MRF.n_t2 + MRF.n_adiasl
    error('check MRF encoding list');
end

if MRF.n_sat > 0
    if MRF.n_sat ~= numel(SAT.sat_rec_time)
        error('number of Saturations is incorrect');
    end
end

if MRF.n_inversion > 0
    if MRF.n_inversion ~= numel(INV.inv_rec_time)
        error('number of Inversions is incorrect');
    end
end

if MRF.n_sl > 0
    if MRF.n_sl ~= numel(SL.tSL)
        error('number of Spin-Lock Preparations is incorrect');
    end
end

if MRF.n_mlev > 0
    if MRF.n_mlev ~= numel(MLEV.t2p_prep_times)
        error('number of MLEV Preparations is incorrect');
    end
end

if MRF.n_t2 > 0
    if MRF.n_t2 ~= numel(T2.prep_times)
        error('number of T2 Preparations is incorrect');
    end
end

if MRF.n_adiasl > 0
    if MRF.n_adiasl ~= numel(ADIASL.adia_prep_times)
        error('number of adiabatic SL Preparations is incorrect');
    end
end


