%% init pulseq
% basic PRESS (Point-RESolved-Spectroscopy) readout
clear
seq_name = 'press';

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
FOV.Nx       = 1;
FOV.Ny       = 1;
FOV.Nz       = 1;
FOV.fov_x    = 10 *1e-3;
FOV.fov_y    = 10 *1e-3;
FOV.fov_z    = 10 *1e-3;
FOV.dx       = FOV.fov_x;
FOV.dy       = FOV.fov_y;
FOV.dz       = FOV.fov_z;
FOV.x_offset = 0 *1e-3;
FOV.y_offset = 0 *1e-3;
FOV.z_offset = 0 *1e-3;
FOV_init();

%% PRESS parameters

% excitation pulses
PRESS.exc.mode = 'sinc';
PRESS.exc.tbw  = 4;
PRESS.exc.t90  = 3 *1e-3;
PRESS.exc.t180 = 3 *1e-3;
PRESS.exc.phz  = pi/2;
PRESS.exc.phy  = PRESS.exc.phz - pi/2;
PRESS.exc.phx  = PRESS.exc.phy + pi;

% crusher gradients
PRESS.crush.n_twists = 8;
PRESS.crush.duration = 2 *1e-3;

% acquisition
PRESS.acq.adcTimeDesired = 10 *1e-3; % [ms]
PRESS.acq.adcBWDesired   = 100 *1e3; % [Hz]

% repetitions
PRESS.NR     = 5;
PRESS.ndummy = 0;
PRESS.Trec   = 1000 *1e-3;

% effective echo time
PRESS.TE.echo_filling_delay = 0 *1e-3; % [s], 0 -> automatic minimization of effective echo time
PRESS.TE.echo_pos_delay     = 0 *1e-3; % [s]

PRESS = PRESS_init(PRESS, FOV, system, 1);

%% magnetization reset
SAT.mode = 'on';
SAT = SAT_init(SAT, FOV, system);

%% create sequence
for loop_NR = 1-PRESS.ndummy : PRESS.NR
    PRESS_add();
    SAT_add();
    seq.addBlock(mr.makeDelay(PRESS.Trec));
end

%% plot sequence diagram
seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();