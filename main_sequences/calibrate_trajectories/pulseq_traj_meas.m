%% init pulseq
% basis: SPI readout
% use for: measurement of gradient trajectories
clear
seq_name = 'traj_';

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

%% import SPI object for trajectory measurement

% load any backup file containing an SPI object
load('Q:/data/Pulseq/Pulseq_Workspace/mgram/260311/260311_0123/backup_260311_0123_workspace.mat')

%% store the original workspace inside the new workspace
PULSEQ_SPI = PULSEQ;
clear PULSEQ;
PULSEQ.PULSEQ_SPI = PULSEQ_SPI;
clear PULSEQ_SPI;
seq_name = [seq_name PULSEQ.PULSEQ_SPI.seq_id];

%% sequence objects for trajectory measurement

% select method
% duyn:    10.1006/jmre.1998.1396
% robison: 10.1002/mrm.27583
TRAJ.method = 'robison';

% number of repetitions for averaging
TRAJ.Nav = 5;

% number of dummy scans
TRAJ.Ndummy = 50;

% recovery time between repetitions
TRAJ.Trec = 150 *1e-3; % [s]

% analyze phase evolution in a thin slice far from the iso-center
TRAJ.slice_thickness = 2 *1e-3;  % [m] slice thickness
TRAJ.slice_offset    = 50 *1e-3; % [m] slice offset

% flip angle
TRAJ.exc_fa = 20 *pi/180;

% reduce stimulation
TRAJ.lim_slew = 0.5;

% init seq objects for trajectory scans
[TRAJ, FOV] = TRAJ_init(TRAJ, PULSEQ.PULSEQ_SPI.SPI, system);

%% add seq ojects

if strcmp(TRAJ.method, 'duyn')
    for loop_xy = [1 3]
        ndummy = TRAJ.Ndummy;
        for loop_NR = 1 : TRAJ.NR            
            loop_traj = loop_xy;
            for loop_av = 1 - ndummy : TRAJ.Nav
                TRAJ_add();
            end
            loop_traj = loop_xy + 1;
            for loop_av = 1 : TRAJ.Nav
                TRAJ_add();
            end
            ndummy = 0;
        end
    end
end

if strcmp(TRAJ.method, 'robison')
    for loop_xy = [1 5]
        ndummy = TRAJ.Ndummy;
        for loop_NR = 1 : TRAJ.NR            
            loop_traj = loop_xy;
            for loop_av = 1 - ndummy : TRAJ.Nav
                loop_traj = loop_xy;
                TRAJ_add();
                loop_traj = loop_xy + 2;
                TRAJ_add();
                loop_traj = loop_xy + 1;
                TRAJ_add();
                loop_traj = loop_xy + 3;
                TRAJ_add();
            end           
            ndummy = 0;
        end
    end
end

%% plot sequence diagram
seq.plot('TimeRange', mr.calcDuration(TRAJ.Trec)*([TRAJ.Ndummy*4, TRAJ.Ndummy*4+10]));

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();