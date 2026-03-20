%% init pulseq
% basic 3D BSSFP readout
clear
seq_name = 'bssfp';

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
FOV.Nx    = 256; % matrix sixe: x
FOV.Ny    = 256; % matrix sixe: y
FOV.Nz    =  48; % matrix sixe: z
FOV.fov_x = 256 *1e-3; % [m] FOV: x direction
FOV.fov_y = 256 *1e-3; % [m] FOV: y direction
FOV.fov_z =  96 *1e-3; % [m] FOV: z direction
FOV.dz    = FOV.fov_z * 0.8; % [m] slab thickness (reduce vs. FOV.fov_z to prevent aliasing)
FOV_init();

%% BSSFP sequence parameters
BSSFP.TR        = 0  *1e-3;   % [s] repetition time, use 0 for auto-minimization
BSSFP.t_acq     = 0.00256;    % [s] acquisition time
BSSFP.t_pre_rew = 2.0 *1e-3;  % [s] duration of prephaser and rewinder gradients
BSSFP.os_mode   = 0;          % [0/1] read oversampling: 0 off, 1 on

% slice excitation
BSSFP.exc_mode       = 'sinc';      % 'sinc' or 'sigpy_SLR'
BSSFP.exc_shape      = 'ex';        % 'st' or 'ex' only for sigpy_SLR case 
BSSFP.exc_time       = 2.0 *1e-3;   % [s] excitation time
BSSFP.exc_tbw        = 4;           % [ ] time bandwidth product 
BSSFP.exc_flipangle  = 60 *pi/180;  % [rad] flip angles

% limit gradients
BSSFP.max_grad = 1/sqrt(3);
BSSFP.max_slew = 1/sqrt(3);

% dummy scans
BSSFP.Ndummy = 250;

[BSSFP, ktraj_adc, ktraj_full] = BSSFP_init(BSSFP, FOV, system);

%% add sequence blocks
seq = mr.Sequence(system);

for loop_dummy = 1:BSSFP.Ndummy
    loop_y = 1;
    loop_z = 1;
    BSSFP_add();
end
loop_dummy = 0;

for loop_z = 1:FOV.Nz
    for loop_y = 1:FOV.Ny
        BSSFP_add();
    end
end

%% plot sequence diagram
seq.plot('TimeRange', (BSSFP.Ndummy - 2 + [0 20]) * BSSFP.TR)

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();