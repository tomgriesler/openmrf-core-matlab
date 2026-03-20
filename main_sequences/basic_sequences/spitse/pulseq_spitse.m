%% init pulseq
% basisc SPI-TSE (spiral-turbo-spin-echo) readout
% -> doi.org/10.1002/mrm.29224
clear
seq_name = 'spitse';

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
FOV.Nxy      = 256;
FOV.fov_xy   = 256 *1e-3;
FOV.dz       = 5 *1e-3;
FOV.z_offset = 0 *1e-3;
FOV_init();

%% spiral sequence parameters
SPITSE.scanmode   = 'run';
SPITSE.segmode    = 'fix';
SPITSE.spmode     = 'cont';
SPITSE.initmode   = 'no';
SPITSE.accmode    = 'vd';
SPITSE.fatsat     = 'off';
SPITSE.seqvar_mod = 'none';
SPITSE.T2prep     = 'off';
SPITSE.EncMode    = 'pd_linear';
SPITSE.plotflag   = '111';
SPITSE.dispflag   = 0;

SPITSE.Trec       = 1000 *1e-3;
SPITSE.TE         = 12 *1e-3;
SPITSE.NEcho      = 2;
SPITSE.tEX        = 2.5 *1e-3;
SPITSE.tRef       = 2.5 *1e-3;
SPITSE.flipref    = 60;
SPITSE.flipflag   = 2;
SPITSE.rfex_phase = 0;
SPITSE.tbw        = 3;
SPITSE.slewfac    = 1/sqrt(3);
SPITSE.gradfac    = 1/sqrt(3);

SPITSE.Ndummy = 1;
SPITSE.Nrep   = 5;

[SPITSE, ktraj_adc, ktraj_full, ktraj_reco] = SPITSE_init(SPITSE, FOV, system);

%% fat saturation and reset
FAT.mode = 'on';
FAT = FAT_init(FAT, FOV, system);

SAT.mode = 'on';
SAT = SAT_init(SAT, FOV, system);

%% create sequence
for loop_rep = 1-SPITSE.Ndummy : SPITSE.Nrep

    % fat saturation
    FAT_add();

    % spiral tse readouts
    [seq] = SPITSE_add(seq, system, FOV, SPITSE, loop_rep);

    % saturation
    SAT_add();

    % recover time
    seq.addBlock(mr.makeDelay(SPITSE.Trec));

end


%% plot sequence diagram
seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();