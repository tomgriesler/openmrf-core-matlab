%% init pulseq
% basis: SPI readout
% use for: T1rho, T2rho or T2 mapping
clear
seq_name = 't1p_mapping';

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

% store spin-lock amplitude to .seq name
fSL = 300;
seq_name = [seq_name '_' num2str(fSL) 'Hz'];

%% FOV geometry
FOV.Nxy      = 192;         % [ ] matrix size
FOV.fov_xy   = 300  *1e-3;  % [m] FOV geometry
FOV.dz       = 5   *1e-3;   % [m] slice thickness
FOV.fov_z    = FOV.dz;
FOV_init();

%% spiral sequence parameters
% params: basic
SPI.mode_2D_3D = '2D';     % '2D', '3D' or '3D_stacked'
SPI.TR         = 0 *1e-3;  % [s]  repetition time, 0 for minimization
SPI.TE         = 0 *1e-3;  % [s]  echo time, 0 for minimization
SPI.NR         = 8;        % [ ]  number of repetitions
SPI.Ndummy     = 0;        % [ ]  initial dummy loops
SPI.adcBW      = 400 *1e3; % [Hz] desired receiver bandwidth

% params: slice excitation
SPI.exc_mode      = 'sinc';      % 'sinc' or 'sigpy_SLR'
SPI.exc_shape     = 'ex';        % only for sigpy: 'st' or 'ex' 
SPI.exc_time      = 2.0 *1e-3;   % [s] excitation time
SPI.exc_tbw       = 4;           % [ ] time bandwidth product
SPI.exc_fa_mode   = 'equal';     % 'equal',  'ramped',  'mrf'  
SPI.exc_flipangle = 90 *pi/180;  % [rad] const FA or start of FA ramp -> set [] for auto mode
SPI.reph_duration = 1.0 *1e-3;   % [s] slice rephaser duration

% params: gradient spoiling & rf spoiling
SPI.spoil_nTwist   = 8;           % [ ] number of 2pi twist in slice
SPI.spoil_rf_mode  = 'quad';      % 'lin' or 'quad'
SPI.spoil_rf_inc   = 117 *pi/180; % [rad] rf phase increment
SPI.spoil_duration = 3.5 *1e-3;   % [s] slice spoiler duration

% k-space geometry params
SPI.geo.design        = 'spiral';
SPI.geo.design_fun    = 'hargreaves';
SPI.geo.N_interleaves = 7.9;
SPI.geo.FOV_coeff     = [1 0];
SPI.geo.Ns            = 1e5;
SPI.geo.ds            = 1e-3;
SPI.geo.flag_rv       = 0;

% gradient & slew rate limitation factors
SPI.lim_exc_slew   = 0.8; % slice excitation
SPI.lim_reph_grad  = 0.9; % slice rephasing
SPI.lim_reph_slew  = 0.6; % slice rephasing
SPI.lim_read_grad  = 0.9; % readout
SPI.lim_read_slew  = 0.8; % readout
SPI.lim_spoil_grad = 0.9; % slice spoiling
SPI.lim_spoil_slew = 0.7; % slice spoiling

% k-space projection params
SPI.proj.mode = 'Equal2Pi'; % 'Equal2Pi' 'RoundGoldenAngle'

% calculate SPI pulse objects
[SPI, ktraj_ref] = SPI_init(SPI, FOV, system, 1);

SPI.Trec = mr.makeDelay(5);

%% params: Spin-Lock

% spin-lock pulses
SL.relax_type = {'T1p'};                  % T1p or T2p or T2
SL.seq_type   = {'BSL'};                  % BSL or CSL or RESL
SL.tSL        = [8 : 8 : 80] *1e-3;       % [s]  SL time
SL.fSL        = fSL * ones(size(SL.tSL)); % [Hz] SL amplitude

% excitation pulses
SL.exc_mode  = 'adiabatic_AHP';  % 'adiabatic_AHP', 'sinc', 'sigpy_SLR' or 'bp'
SL.exc_time  = 3.0 *1e-3;        % [s] excitation time
SL.adia_wmax = 600 * 2*pi;       % [rad/s] amplitude of adiabatic pulse

% refocusing pulses
SL.rfc_mode = 'bp';              % 'bp', 'sinc', 'sigpy_SLR' or 'comp'
SL.rfc_time = 1.0 *1e-3;         % [s] refocusing time

SL = SL_init(SL, FOV, system);

%% fat saturation and reset
FAT.mode = 'on';
FAT = FAT_init(FAT, FOV, system);

SAT.mode = 'on';
SAT = SAT_init(SAT, FOV, system);

%% create sequence

ndummy = SPI.Ndummy;

for loop_SL = 1 : SL.nSL    
for loop_NR = 1-ndummy : SPI.NR

    ndummy = 0;

    % saturtion or crusher
    SAT_add();   
    
    % recovery time
    seq.addBlock(SPI.Trec);

    % fat saturation
    FAT_add();

    % sl preparation
    SL_add();

    % spiral imaging
    SPI_add();  

end    
end

%% plot sequence diagram
seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();