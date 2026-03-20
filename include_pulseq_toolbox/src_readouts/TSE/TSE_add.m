% ---------------------------------------------------------
% ------------------ add pulseq objects -------------------
% ------------ readout: Turbo-Spin-Echo (TSE) -------------
% ---------------------------------------------------------

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% ----- sequence loop counter: -----
% loop_TE -> different echos
% loop_TR -> repetitions

%% add sequence blocks
TSE.rf_rfc.phaseOffset = TSE.rfc_phase;
seq.addBlock(TSE.GS1);
seq.addBlock(TSE.GS2, TSE.rf_exc);
seq.addBlock(TSE.GS3, TSE.GR3);
for loop_TE = 1 : TSE.n_echo
    seq.addBlock(TSE.GS4, TSE.rf_rfc);
    TSE.rf_rfc.phaseOffset = TSE.rf_rfc.phaseOffset + pi;    
    if (loop_TR>0)
        seq.addBlock(TSE.GS5, TSE.GR5, TSE.GPpre(loop_TR,loop_TE).GPpre);
        seq.addBlock(TSE.GR6, TSE.adc);
        seq.addBlock(TSE.GS7, TSE.GR7, TSE.GPrew(loop_TR,loop_TE).GPrew);
    else
        seq.addBlock(TSE.GS5, TSE.GR5, TSE.GPpre(1,1).GPpre);
        seq.addBlock(TSE.GR6);
        seq.addBlock(TSE.GS7, TSE.GR7, TSE.GPrew(1,1).GPrew);
    end    
end
seq.addBlock(TSE.GS4);
seq.addBlock(TSE.GS5);
seq.addBlock(TSE.tr_filling_delay);
