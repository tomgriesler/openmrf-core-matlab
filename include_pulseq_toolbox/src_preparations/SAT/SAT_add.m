% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% add magnetization reset, saturation pulse + crusher
if ~strcmp(SAT.mode, 'off')
    if ~exist('loop_SAT', 'var')
        loop_SAT = 1;
    end   
    seq.addBlock(SAT.gx_crush);
    seq.addBlock(SAT.rf);
    seq.addBlock(SAT.gy_crush, SAT.gz_crush);
    seq.addBlock(SAT.sat_rec_delay(loop_SAT));
end