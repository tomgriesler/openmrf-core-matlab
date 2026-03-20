function params = setupParameters_LowrankMRF2D( )

    params = struct();
    params.numIter          = 50;     % total number of iterations
    params.block_dim        = [6,6];  % locally low-rank patch size: should divide evenly into the matrix size (e.g., 192 is divisible by 6)
    params.block_step       = 6;      % overlap between local patches (if equal to patch size above, then patches are non-overlapping)
    params.lambdaLLR        = 0.02;   % locally low-rank regularization
    params.lambdaSpatialTV  = 0.003;  % total variation regularization
    params.lambdaWav        = 0;      % wavelet regularization
    params.lambdaTemporalTV = 0;      % total variation denoising: MRF contrast dimension
    params.lambdaTemporalFT = 0;      % total variation denoising: time dimension
    params.betaMethod       = 'Dai-Yuan';
    params.beta             = 0.6;
    params.alpha            = 0.01;   % adjusts the OGM1 step size
    params.stopThresh       = 0;      % stopping criterion
    params.updateFreq       = 0;      % how often to display intermediate images  
    params.t0               = 1;

end