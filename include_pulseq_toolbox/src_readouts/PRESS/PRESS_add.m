% ---------------------------------------------------------
% ------------------ add pulseq objects -------------------
% ----- readout: Point RESolved Spectroscopy (PRESS) ------
% ---------------------------------------------------------

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% ----- sequence loop counter: -----
% loop_NR -> repetitions

%% add sequence blocks

seq.addBlock(PRESS.exc.rf90z, PRESS.exc.g90z);
seq.addBlock(PRESS.exc.g90z_reph);
seq.addBlock(PRESS.TE.delay1);

seq.addBlock(PRESS.crush.gy);
seq.addBlock(PRESS.exc.rf180y, PRESS.exc.g180y);
seq.addBlock(PRESS.crush.gy);
seq.addBlock(PRESS.TE.delay2);

seq.addBlock(PRESS.crush.gx);
seq.addBlock(PRESS.exc.rf180x, PRESS.exc.g180x);
seq.addBlock(PRESS.crush.gx);
seq.addBlock(PRESS.TE.delay3);

if loop_NR>0
    seq.addBlock(PRESS.acq.adc, PRESS.acq.delay);
else
    seq.addBlock(PRESS.acq.delay);
end