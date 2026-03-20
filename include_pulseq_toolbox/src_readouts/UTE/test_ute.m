clear

pulseq_init();

FOV.fov_xy   = 256 *1e-3;
FOV.Nxy      = 256;
FOV.dz       = 5 *1e-3;
FOV.z_offset = 0;
FOV_init();

% excitation pulse
UTE.exc_fa    = 10 *pi/180;
UTE.exc_time  = 0.5 *1e-3;
UTE.exc_tbw   = 2;
UTE.exc_phase = pi/2;
UTE.exc_mode  = 'sinc';

% porjections
UTE.NR       = 144;         % numer of repetitions / spokes
UTE.TR       = 20 *1e-3;    % [s] repetition time
UTE.TE       = 70 *1e-6;    % [s] echo time
UTE.Ndummy   = 10;          % number of dummy scans
UTE.phi_mode = 'Equal2Pi';  % 'EqualPi',  'GoldenAngle',  'RandomPi',  'RandomGoldenAngle'

% adc
UTE.adcTime    = 0.52 *1e-3; % [s] adc duration
UTE.ro_os      = 1;          % readout oversampling
UTE.ro_discard = 0;          % dummy ADC samples to discard (due to ADC filter)

% gradient spoiling
UTE.spoil_ro = 1;        % additional k-max excursion for readout spoiling

% rf spoiling
UTE.spoil_rf_mode = 'quad';       % rf spoiling mode: 'lin' or 'rad'
UTE.spoil_rf_inc  = 117 *pi/180;  % rf spoiling increment [rad]

[UTE, ktraj_adc, ktraj_full, ktraj_reco] = UTE_init(UTE, FOV, system);

%% start the sequence
for loop_NR = (1-UTE.Ndummy) : UTE.NR
    UTE_add();
end

seq.plot()