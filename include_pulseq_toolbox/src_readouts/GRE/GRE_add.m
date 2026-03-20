% ---------------------------------------------------------
% ------------------ add pulseq objects -------------------
% ------------- readout: GRadient-Echo (GRE) --------------
% ---------------------------------------------------------

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% ----- sequence loop counter: -----
% loop_Ny -> phase encoding steps
% loop_rf_inc -> rf spoiling
% loop_FA -> different flip angles
% loop_TE -> different echo times

%% init loop counters
if ~exist('loop_rf_inc', 'var')
    loop_rf_inc = 0;
end
if ~exist('loop_FA', 'var')
    loop_FA = 1;
end
if ~exist('loop_TE', 'var')
    loop_TE = 1;
end
if loop_Ny<1
    temp_loop = 1;
else
    temp_loop = loop_Ny;
end

%% rf spoilng
temp_phase                  = 0 + GRE.spoil_rf_inc * loop_rf_inc ^ GRE.spoil_rf_pow;
GRE.rf(loop_FA).phaseOffset = mod(temp_phase,        2*pi);
GRE.adc.phaseOffset         = mod(temp_phase + pi/2, 2*pi);
loop_rf_inc                 = loop_rf_inc + 1;

%% add sequence blocks
seq.addBlock(GRE.rf(loop_FA), GRE.gz);
seq.addBlock(GRE.gx_pre, GRE.gy_pre(temp_loop), GRE.gz_reph);
seq.addBlock(mr.makeDelay(GRE.delayTE(loop_TE)));
if loop_Ny<1
    seq.addBlock(GRE.gx); % dummy loop
else
    seq.addBlock(GRE.gx, GRE.adc);
end
seq.addBlock(GRE.gx_rew, GRE.gy_rew(temp_loop), GRE.gz_spoil);
seq.addBlock(mr.makeDelay(GRE.delayTR(loop_TE)))

clear temp_loop temp_phase;
