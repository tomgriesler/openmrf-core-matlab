%% init pulseq
% basic 3D cones readout
clear
seq_name = '3D_cones';

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
FOV.Nxy    = 256;             % [ ] matrix size
FOV.fov_xy = 256  *1e-3;      % [m] FOV geometry
FOV.fov_z  = FOV.fov_xy;      % [m] FOV geometry; 3D readout is isotropic
FOV.dz     = FOV.fov_z * 0.8; % [m] slab thickness (reduce to prevent aliasing)
FOV_init();

%% spiral sequence parameters

% params: basic
SPI.mode_2D_3D = '3D';     % '2D', '3D' or '3D_stacked'
SPI.TR         = 25 *1e-3; % [s]  repetition time, 0 for minimization
SPI.TE         = 0 *1e-3;  % [s]  echo time, 0 for minimization
SPI.NR         = 10000;    % [ ]  number of repetitions
SPI.Ndummy     = 250;      % [ ]  initial dummy loops
SPI.adcBW      = 400 *1e3; % [Hz] desired receiver bandwidth

% params: slice excitation
SPI.exc_mode      = 'sinc';     % 'sinc' or 'sigpy_SLR'
SPI.exc_shape     = 'ex';       % only for sigpy: 'st' or 'ex' 
SPI.exc_time      = 1.0 *1e-3;  % [s] excitation time
SPI.exc_tbw       = 6;          % [ ] time bandwidth product
SPI.exc_fa_mode   = 'equal';    % 'equal',  'ramped',  'mrf'  
SPI.exc_flipangle = 10 *pi/180; % [rad] const FA or start of FA ramp -> set [] for auto mode
SPI.reph_duration = 0.5 *1e-3;  % [s] slice rephaser duration

% params: gradient spoiling & rf spoiling
SPI.spoil_nTwist   = 2 * FOV.Nxy; % [ ] number of 2pi twist in slice
SPI.spoil_rf_mode  = 'quad';      % 'lin' or 'quad'
SPI.spoil_rf_inc   = 117 *pi/180; % [rad] rf phase increment
SPI.spoil_duration = 5 *1e-3;     % [s] duration of slice spoiler

% k-space geometry params
SPI.geo.design    = 'cones';
SPI.geo.N_loops   = 16;
SPI.geo.theta     = pi/8;
SPI.geo.FOV_coeff = [1 -0.5];
SPI.geo.Ns        = 1e5;
SPI.geo.ds        = 1e-3;
SPI.geo.flag_rv   = 0;

% gradient & slew rate limitation factors
SPI.lim_exc_slew   = 0.8; % slice excitation
SPI.lim_reph_grad  = 0.9; % slice rephasing
SPI.lim_reph_slew  = 0.4; % slice rephasing
SPI.lim_read_grad  = 0.9; % readout
SPI.lim_read_slew  = 0.5; % readout
SPI.lim_spoil_grad = 0.8; % slice spoiling
SPI.lim_spoil_slew = 0.6; % slice spoiling

% k-space projection params
SPI.proj.mode = 'FibonacciSphere'; % 'Equal2Pi' 'RoundGoldenAngle'

% calculate SPI pulse objects
[SPI, ktraj_ref] = SPI_init(SPI, FOV, system, 1);

%% noise pre-scans
SPI.Nnoise = 16;
SPI_add_prescans();

%% create sequence
for loop_NR = 1-SPI.Ndummy : SPI.NR
    SPI_add();
end

%% plot sequence diagram
% seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();