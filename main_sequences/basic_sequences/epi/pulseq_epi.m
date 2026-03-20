%% init pulseq
% basic EPI (echo-planar-imaging) readout
clear
seq_name = 'epi';

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
FOV.Nxy      = 64;
FOV.fov_xy   = 256 *1e-3;
FOV.dz       = 5 *1e-3;
FOV.z_offset = 0 *1e-3;
FOV_init();

%% EPI sequence paramsters
EPI.pe_enable         = 1;   % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
EPI.ro_os             = 1;   % oversampling factor (in contrast to the product sequence we don't really need it)
EPI.readoutTime       = 4.2 *1e-4;   % this controls the readout bandwidth
EPI.partFourierFactor = 1;   % partial Fourier factor: 1: full sampling 0: start with ky=0

EPI.exc_time = 3 *1e-3;
EPI.exc_tbw  = 4;
EPI.exc_mode = 'sinc';

[EPI, ktraj_adc, ktraj_full, ktraj_reco] = EPI_init(EPI, FOV, system);

%% add sequence blocks
for j=1:10
    EPI_add();
    seq.addBlock(mr.makeDelay(1));
end

%% plot sequence diagram
seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();