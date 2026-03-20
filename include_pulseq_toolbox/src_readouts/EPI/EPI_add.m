% ---------------------------------------------------------
% ------------------ add pulseq objects -------------------
% -------------- readout: echo-planar (EPI) ---------------
% ---------------------------------------------------------

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% ----- sequence loop counter: -----
% loop_Ny -> phase encoding steps

%% add sequence blocks
temp_gx = EPI.gx;
seq.addBlock(EPI.rf, EPI.gz);
seq.addBlock(EPI.gxPre, EPI.gyPre, EPI.gzReph);
for loop_Ny = 1 : EPI.Ny_meas
    if loop_Ny==1
        seq.addBlock(temp_gx, EPI.gy_blipup, EPI.adc); % Read the first line of k-space with a single half-blip at the end
    elseif loop_Ny == EPI.Ny_meas
        seq.addBlock(temp_gx, EPI.gy_blipdown, EPI.adc); % Read the last line of k-space with a single half-blip at the beginning
    else
        seq.addBlock(temp_gx, EPI.gy_blipdownup, EPI.adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
    end
    temp_gx.amplitude = -temp_gx.amplitude;   % Reverse polarity of read gradient
end
clear temp_gx;