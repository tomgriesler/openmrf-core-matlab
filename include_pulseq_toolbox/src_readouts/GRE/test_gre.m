clear

pulseq_init();

%% FOV geometry
FOV.Nx       = 128;
FOV.Ny       = 128;
FOV.fov_x    = 240 *1e-3;
FOV.fov_y    = 240 *1e-3;
FOV.dz       = 5 *1e-3;
FOV.z_offset = 0 *1e-3;
FOV_init();

%% GRE sequence parameters
GRE.TEs      = [2 4 6 8] *1e-3;  % [s] echo times, use 0 for auto-minimization
GRE.TRs      = [10] *1e-3;     % [s] repetition time, use 0 for auto-minimization
GRE.t_acq    = 3.2 *1e-3;    % [s] acquisition time
GRE.os_mode  = 0;            % read oversampling: 0 off, 1 on
GRE.t_pre    = 1.5 *1e-3;    % [s] time for prehaser and rephaser gradients, use [] for auto-minimization
GRE.t_spoil  = 1.5 *1e-3;    % [s] time for rewinder and spoiler  gradients, use [] for auto-minimization
GRE.n_rep    = 1;            % [ ] number of repetitions, averages
GRE.ndummy   = 50;           % [ ] number of dummy scans

% slice excitation
GRE.exc_mode       = 'sinc';        % 'sinc' or 'sigpy_SLR'
GRE.exc_shape      = 'ex';          % 'st' or 'ex' only for sigpy_SLR case 
GRE.exc_time       = 2.0 *1e-3;     % [s] excitation time
GRE.exc_tbw        = 6;             % [ ] time bandwidth product 
GRE.exc_flipangle  = [20] *pi/180;  % [rad] flip angles

% slice spoiling
GRE.spoil_nTwist  = 8;  % [ ] number of 2pi twists in z-direction, 0 for balanced moments

% spoiling
GRE.spoil_rf_mode = 'quad';        % lin or quad
GRE.spoil_rf_inc  = 117 *pi/180;   % [deg] use 117 for flash-type and 180 for ssfp-type

[GRE, ktraj_adc, ktraj_full] = GRE_init(GRE, FOV, system);