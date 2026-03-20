% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% CEST preparation
if (loop_ppm < CEST.n_ppm + 1)
    CEST.rf.freqOffset  = CEST.f_off(loop_ppm);
    seq.addBlock(CEST.rf);
else
    seq.addBlock(mr.makeDelay(CEST.ref_delay));
end
seq.addBlock(CEST.gx_crush, CEST.gy_crush, CEST.gz_crush);