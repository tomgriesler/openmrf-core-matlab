%% init pulseq
% basic GRE (gradient-echo) readout
clear
seq_name = 'gre';

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
FOV.Nx       = 256;
FOV.Ny       = 256;
FOV.fov_x    = 256 *1e-3;
FOV.fov_y    = 256 *1e-3;
FOV.dz       = 5 *1e-3;
FOV.z_offset = 0 *1e-3;
FOV_init();

%% GRE sequence parameters
GRE.TEs      = [4] *1e-3;    % [s] echo times
GRE.TRs      = 50 *1e-3;     % [s] repetition time, use 0 for auto-minimization
GRE.t_acq    = 3.2 *1e-3;    % [s] acquisition time
GRE.os_mode  = 0;            % read oversampling: 0 off, 1 on
GRE.t_pre    = 1.5 *1e-3;    % [s] time for prehaser and rephaser gradients, use [] for auto-minimization
GRE.t_spoil  = 1.5 *1e-3;    % [s] time for rewinder and spoiler  gradients, use [] for auto-minimization
GRE.n_rep    = 1;            % [ ] number of repetitions, averages
GRE.ndummy   = 50;           % [ ] number of dummy scans

% slice excitation
GRE.exc_mode       = 'sinc';        % 'sinc' or 'sigpy_SLR'
GRE.exc_shape      = 'ex';          % 'st' or 'ex' only for sigpy_SLR case 
GRE.exc_time       = 1.0 *1e-3;     % [s] excitation time
GRE.exc_tbw        = 4;             % [ ] time bandwidth product 
GRE.exc_flipangle  = [20] *pi/180;  % [rad] flip angles

% spoiling
GRE.spoil_rf_mode = 'quad';        % lin or quad
GRE.spoil_rf_inc  = 117 *pi/180;   % [deg] use 117 for flash-type and 180 for ssfp-type
GRE.spoil_nTwist  = 8;             % [ ] number of 2pi twists in z-direction, 0 for balanced

[GRE, ktraj_adc, ktraj_full] = GRE_init(GRE, FOV, system);

%% add sequence blocks
for loop_FA = 1 : GRE.n_FAs
for loop_TE = 1 : GRE.n_TEs

    ndummy = GRE.ndummy;
    for loop_rep = 1:GRE.n_rep
        for loop_Ny  = -ndummy+1 : FOV.Ny    
            GRE_add();
        end
        ndummy = 0;
    end

end 
end

%% plot sequence diagram
seq.plot('TimeRange',[GRE.ndummy+1 GRE.ndummy+5]*GRE.TRs(1))

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();