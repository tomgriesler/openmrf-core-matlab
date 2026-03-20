% ---------------------------------------------------------
% ------------------ add pulseq objects -------------------
% ------------ readout: Ultra-short TE (UTE) --------------
% ---------------------------------------------------------

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% ----- loop counter: -----
% loop_rf_inc -> rf spoiling
% loop_NR -> repetitions

%% spoiling and interleaving
if ~exist('loop_rf_inc', 'var')
    loop_rf_inc = 0;
end
if strcmp(UTE.spoil_rf_mode, 'lin')
    temp_phase = UTE.exc_phase + UTE.spoil_rf_inc * loop_rf_inc;
elseif strcmp(UTE.spoil_rf_mode, 'quad')
    temp_phase = UTE.exc_phase + UTE.spoil_rf_inc * loop_rf_inc^2;
else
    error('spoiling increment mode unknown!');
end
UTE.rf.phaseOffset  = mod(temp_phase,        2*pi);
UTE.adc.phaseOffset = mod(temp_phase - pi/2, 2*pi);
temp_phi            = UTE.phi(loop_NR + UTE.Ndummy);
loop_rf_inc         = loop_rf_inc + 1;

%% add sequence blocks
% 1st readout: negative slice gradient
if (loop_NR>0)
    seq.addBlock(mr.rotate('z', temp_phi, UTE.rf, UTE.gz_neg, UTE.gx, UTE.adc));
else
    seq.addBlock(mr.rotate('z', temp_phi, UTE.rf, UTE.gz_neg, UTE.gx));
end
seq.addBlock(UTE.delayTR);

% spoiling and interleaving
if strcmp(UTE.spoil_rf_mode, 'lin')
    temp_phase = UTE.exc_phase + UTE.spoil_rf_inc * loop_rf_inc;
elseif strcmp(UTE.spoil_rf_mode, 'quad')
    temp_phase = UTE.exc_phase + UTE.spoil_rf_inc * loop_rf_inc^2;
else
    error('spoiling increment mode unknown!');
end
UTE.rf.phaseOffset  = temp_phase;
UTE.adc.phaseOffset = temp_phase - pi/2;
temp_phi            = UTE.phi(loop_NR + UTE.Ndummy);
loop_rf_inc         = loop_rf_inc + 1;

% 2nd readout: positive slice gradient
if (loop_NR>0)
    seq.addBlock(mr.rotate('z', temp_phi, UTE.rf, UTE.gz_pos, UTE.gx, UTE.adc));
else
    seq.addBlock(mr.rotate('z', temp_phi, UTE.rf, UTE.gz_pos, UTE.gx));
end
seq.addBlock(UTE.delayTR);
clear temp_phase temp_phi;