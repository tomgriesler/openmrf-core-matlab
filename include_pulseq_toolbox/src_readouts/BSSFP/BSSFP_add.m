% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% ---------------------------------------------------------
% ------------------ add pulseq objects -------------------
% ------------------ readout: 3D BSSFP --------------------
% ---------------------------------------------------------

% ----- sequence loop counter: -----
% loop_y      -> phase encoding steps
% loop_z      -> partition encoding steps
% loop_FA     -> different flip angles
% loop_rf_inc -> rf spoiling

%% init loop counters
if ~exist('loop_rf_inc', 'var')
    loop_rf_inc = 0;
end
if ~exist('loop_FA', 'var')
    loop_FA = 1;
end
if ~exist('loop_dummy', 'var')
    loop_dummy = 0;
end

%% rf spoilng
temp_phase                    = pi * loop_rf_inc;
BSSFP.rf(loop_FA).phaseOffset = mod(temp_phase,        2*pi);
BSSFP.adc.phaseOffset         = mod(temp_phase + pi/2, 2*pi);
loop_rf_inc                   = loop_rf_inc + 1;
clear temp_phase;

%% add sequence blocks
seq.addBlock(BSSFP.rf(loop_FA), BSSFP.gz);
seq.addBlock(BSSFP.tr_fill);
seq.addBlock(BSSFP.gx_pre_rew, BSSFP.gy_pre(loop_y), BSSFP.gz_pre(loop_z));
if loop_dummy==0
    seq.addBlock(BSSFP.gx, BSSFP.adc);
else
    seq.addBlock(BSSFP.gx);
end
seq.addBlock(BSSFP.gx_pre_rew, BSSFP.gy_rew(loop_y), BSSFP.gz_rew(loop_z));
seq.addBlock(BSSFP.tr_fill);
