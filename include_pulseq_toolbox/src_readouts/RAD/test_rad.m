clear

pulseq_init();

FOV.fov_xy   = 256e-3;
FOV.Nxy      = 144;
FOV.dz       = 6 *1e-3;
FOV.z_offset = 0;
FOV_init();

% excitation pulse
RAD.exc_fa    = 10 *pi/180;
RAD.exc_time  = 2 *1e-3;
RAD.exc_tbw   = 6;
RAD.exc_phase = pi/2;
RAD.exc_mode  = 'sinc';

% porjections
RAD.NR       = 128;         % numer of repetitions / spokes
RAD.TR       = 10 *1e-3;    % [s] repetition time
RAD.Ndummy   = 10;          % number of dummy scans
RAD.phi_mode = 'Equal2Pi';  % 'Equal2Pi',  'GoldenAngle',  'Random2Pi',  'RandomGoldenAngle'

% adc
RAD.adcTime = 0.6 *1e-3; % [s] adc duration
RAD.ro_os   = 2;         % readout oversampling

% gradient spoiling
RAD.spoil_sl = 8;        % spoil area compared to the slice thickness
RAD.mode_rewind = 'on';  % add rewinder gradient during slice spoiling
RAD.spoil_ro = 4;        % additional k-max excursion for readout spoiling

% rf spoiling
RAD.spoil_rf_mode = 'quad';       % rf spoiling mode: 'lin' or 'rad'
RAD.spoil_rf_inc  = 117 *pi/180;  % rf spoiling increment [rad]

% init RAD objects
[RAD, ktraj_adc, ktraj_full, ktraj_reco] = RAD_init(RAD, FOV, system);

%% start the sequence
for loop_NR = (1-RAD.Ndummy) : RAD.NR
    RAD_add();
end

