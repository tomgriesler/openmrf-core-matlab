%% init pulseq
% basic RAD (radial gradient-echo) readout
clear
seq_name = 'rad_basic';

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
RAD.exc_fa    = 20 *pi/180;
RAD.exc_time  = 1  *1e-3;
RAD.exc_tbw   = 6;
RAD.exc_phase = pi/2;
RAD.exc_mode  = 'sinc';

% porjections
RAD.NR       = 233;         % numer of repetitions / spokes
RAD.TR       = 50 *1e-3;    % [s] repetition time
RAD.Ndummy   = 50;          % number of dummy scans
RAD.phi_mode = 'GoldenAngle';  % 'EqualPi',  'GoldenAngle',  'RandomPi',  'RandomGoldenAngle'

% adc
RAD.adcTime = 0.6 *1e-3; % [s] adc duration
RAD.ro_os   = 2;         % readout oversampling

% gradient spoiling
RAD.spoil_sl = 8;        % spoil area compared to the slice thickness
RAD.mode_rewind = 'on';  % add rewinder gradient during slice spoiling
RAD.spoil_ro = 0;        % additional k-max excursion for readout spoiling

% rf spoiling
RAD.spoil_rf_mode = 'quad';       % rf spoiling mode: 'lin' or 'rad'
RAD.spoil_rf_inc  = 117 *pi/180;  % rf spoiling increment [rad]

% init RAD objects
[RAD, ktraj_adc, ktraj_full, ktraj_reco] = RAD_init(RAD, FOV, system);

%% add sequence blocks
for loop_NR = (1-RAD.Ndummy) : RAD.NR
    RAD_add();
end

%% plot sequence diagram
seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();