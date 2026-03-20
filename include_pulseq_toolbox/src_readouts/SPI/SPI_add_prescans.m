% ---------------------------------------------------------
% ------------------ add pulseq objects -------------------
% ----------------- readout: spiral (SPI) -----------------
% -------------------- Noise pre-scans --------------------
% ---------------------------------------------------------

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

if isfield(SPI, 'Nnoise')
    for loop_noise = 1:SPI.Nnoise
        seq.addTRID('noise_prescan');
        seq.addBlock(SPI.adc, ceil((mr.calcDuration(SPI.adc) + 1e-3) / system.blockDurationRaster) * system.blockDurationRaster);
    end
end
