% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% WASABI preparation
if (loop_ppm < WASABI.n_ppm + 1)
    WASABI.rf.freqOffset  = WASABI.f_off(loop_ppm);
    seq.addBlock(WASABI.rf);
else
    seq.addBlock(mr.makeDelay(WASABI.ref_delay));
end
seq.addBlock(WASABI.gx_crush, WASABI.gy_crush, WASABI.gz_crush);