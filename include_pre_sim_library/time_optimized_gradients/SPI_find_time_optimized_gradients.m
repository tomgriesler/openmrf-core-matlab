function [C, time, g, s, k, phi, sta, stb] = SPI_find_time_optimized_gradients(C, rv, g0, gfin, gmax, smax, dt, ds)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% combine all inputs in a single 1d array
inputs_all = [C(:); rv+0; g0; gfin; gmax; smax; dt; ds];

% create hash
opt_grad_hash = [pulseq_get_path('SPI_find_time_optimized_gradients') 'opt_grad_' pulseq_get_wave_hash(inputs_all) '.mat'];

% find optimized waveform  
if isfile(opt_grad_hash)
    load(opt_grad_hash);
    phi  = [];
    sta  = [];
    stb  = [];
else
    C    = [];
    time = [];
    g    = [];
    s    = [];
    k    = [];
    phi  = [];
    sta  = [];
    stb  = [];
end

end