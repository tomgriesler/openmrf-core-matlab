% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% add magnetization inversion pulse
if ~strcmp(ADIASL.mode, 'off')
    for temp_n_hs = 1:ADIASL.N_HS(loop_ADIASL)
        temp_rf = ADIASL.rf;
        if mod(temp_n_hs,2)==1        
            temp_rf.signal = ADIASL.f1 .* exp(1i*ADIASL.phi_down);
        else
            temp_rf.signal = ADIASL.f1 .* exp(1i*ADIASL.phi_up);
        end
        temp_rf.phaseOffset = (ADIASL.phi + ADIASL.phi_list(temp_n_hs)*pi);        
        seq.addBlock(temp_rf);    
    end
    seq.addBlock(ADIASL.gx_crush, ADIASL.gy_crush, ADIASL.gz_crush);
    clear temp_n_hs temp_rf;
end