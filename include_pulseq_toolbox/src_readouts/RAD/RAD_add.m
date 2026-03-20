% ---------------------------------------------------------
% ------------------ add pulseq objects -------------------
% ----------------- readout: radial (RAD) -----------------
% ---------------------------------------------------------

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% ----- sequence loop counter: -----
% loop_NR -> repetitions
% loop_rf_inc -> rf spoiling

%% rf spoiling
if ~exist('loop_rf_inc', 'var')
    loop_rf_inc = 0;
end
temp_phase          = 0 + RAD.spoil_rf_inc * loop_rf_inc ^RAD.spoil_rf_pow;
RAD.rf.phaseOffset  = mod(temp_phase,        2*pi);
RAD.adc.phaseOffset = mod(temp_phase + pi/2, 2*pi);
loop_rf_inc         = loop_rf_inc + 1;

%% radial interleaving
if (loop_NR>0)
    temp_phi = RAD.phi(loop_NR);
else
    temp_phi = RAD.phi(1);
end

%% add sequence blocks
seq.addBlock(mr.rotate('z', temp_phi, RAD.rf, RAD.gzComb, RAD.gxPre));
if (loop_NR>0)
    seq.addBlock(mr.rotate('z', temp_phi, RAD.gxComb, RAD.gzSpoil, RAD.adc));
else
    seq.addBlock(mr.rotate('z', temp_phi, RAD.gxComb, RAD.gzSpoil));
end
seq.addBlock(RAD.tr_filling_delay);
clear temp_phi temp_phase;