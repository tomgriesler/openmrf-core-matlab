%% init pulseq
% basic UTE (ultra-short-echo_time) readout
clear
seq_name = 'ute';

% main flags
flag_backup = 0; % 0: off,  1: only backup,  2: backup and send .seq
flag_report = 0; % 0: off,  1: only timings, 2: full report (slow)
flag_pns    = 0; % 0: off,  1: simulate PNS stimulation

% optional: select scanner
% pulseq_scanner = 'Siemens_Sola_1,5T_MIITT';

% optional: select pns sim orientation
% pns_orientation = 'coronal';

% init system, seq object and load pulseq user information
pulseq_init();

%% FOV geometry
FOV.fov_xy   = 256 *1e-3;
FOV.Nxy      = 128;
FOV.dz       = 5 *1e-3;
FOV.z_offset = 0;
FOV_init();

%% RAD sequence paramsters

% excitation pulse
UTE.exc_fa    = 20 *pi/180;
UTE.exc_time  = 0.5 *1e-3;
UTE.exc_tbw   = 2;
UTE.exc_phase = pi/2;
UTE.exc_mode  = 'sinc';

% porjections
UTE.NR       = 233;         % numer of repetitions / spokes
UTE.TR       = 50 *1e-3;    % [s] repetition time
UTE.TE       = 70 *1e-6;    % [s] echo time
UTE.Ndummy   = 50;          % number of dummy scans
UTE.phi_mode = 'GoldenAngle';  % 'EqualPi',  'GoldenAngle',  'RandomPi',  'RandomGoldenAngle'

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

%% add sequence blocks
for loop_NR = (1-UTE.Ndummy) : UTE.NR
    UTE_add();
end

%% plot sequence diagram
seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();