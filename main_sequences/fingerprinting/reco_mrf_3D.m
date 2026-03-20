%% load mrf study
clear

% enter paths of mrf and trajectory rawdata
% note: GE, United Imaging and Philips need explicit description of the
% pulseq backup paths. for Siemens, rawdata and pulseq backups can be
% automatically linked

vendor           = 'Siemens';
path_raw_mrf     = 'Q:/data/Pulseq/Rawdata/mgram/SolaMIITT/250311_tests_v151/meas_MID00317_FID84194_260311_0124_mgram_mrf_3D_stackedyun.dat';
path_raw_traj    = 'Q:/data/Pulseq/Rawdata/mgram/SolaMIITT/250311_tests_v151/meas_MID00321_FID84198_260311_0203_mgram_traj_260311_0123.dat';
path_backup_mrf  = [];
path_backup_traj = [];

%% load MRF rawdata
[DATA, NOISE, PULSEQ, study_info] = pulseq_read_meas(path_raw_mrf, path_backup_mrf, vendor);

%% sort kz partitions and use fft in z direction (not tested)

% read dimensions, convert to single, restructure rawdata
[NCoils, Nseg, NRead] = size(DATA);
Nz   = PULSEQ.FOV.Nz;
NR   = PULSEQ.SPI.NR;
Nseg = Nseg / NR;
DATA = single(DATA);
DATA = permute(reshape(permute(DATA, [2, 3, 1]), [NR, Nseg, NRead, NCoils]), [2, 1, 3, 4]);

% sort kz partitions and average repetitions
DATA_3D    = complex(zeros(Nz, NR, NRead, NCoils, 'single'));
kz_missing = [];
for j = 1:Nz
    ind = find(PULSEQ.SPI.kz_table==j);
    if isempty(ind)
        warning(['missing k-space partition no.: ' num2str(j)]);
        kz_missing = [kz_missing; j]; % store missing partitions -> GRAPPA?
    else
        temp_data = DATA(ind, :,:,:);
        if size(temp_data,1) == 1
            temp_data = squeeze(temp_data);
        else
            temp_data = squeeze(mean(temp_data));
        end
        DATA_3D(j,:,:,:) = temp_data;
    end
end
clear j ind temp_data DATA;

% fft along kz direction
DATA_3D_FFTz = kspace2image(DATA_3D, [1, 0, 0, 0]);
DATA_3D_FFTz = permute(DATA_3D_FFTz, [1, 4, 2, 3]);
clear DATA_3D;

%% load measured trajectory
if ~isempty(path_raw_traj)
    [ktraj_meas, ktraj_hash] = TRAJ_reco(path_raw_traj, path_backup_traj, vendor, [8:8:48]);
    if ~strcmp(ktraj_hash, pulseq_get_wave_hash(SPI_load_ktraj(PULSEQ)))
        ktraj_meas = []; % ensure compatibility of trajectory scan
    end
else
    ktraj_meas = [];
end

%% parameters: dictionary simulation
params_dict.seq_path    = [];       % full path of .seq file; [] for automatic search in user backups
params_dict.sim_mode    = 'BLOCH';  % 'BLOCH' or 'EPG'
params_dict.comp_energy = 0.9999;   % svd compression energy: 0 for uncompressed dictionary
params_dict.N_iso       = 1000;     % number of isochromats for bloch simulation
params_dict.s_fac       = 1;        % factor for out-of-slice simulation -> no out-of slice simulation for stacked 3D!
params_dict.f0          = [];       % larmor frequency f0: required to simulate rf pulses with ppm freq offsets
params_dict.time_stamps = [];       % adc times stamps: required for correction of trigger delays in cMRF
params_dict.soft_delays = [];       % soft delay user input: required for correction of the acq window in cMRF
params_dict.flag_kz     = 1;        % find kz partitions for stacked 3D MRF -> eliminate unnecessary partitions
params_dict.echo_mode   = [];       % echo mode; default: 'spiral_out'

%% parameters: dictionary and look-up table

% T1, T2
P.T1.range = [0.01,  4]; P.T1.factor = 1.05;
P.T2.range = [0.001, 3]; P.T2.factor = 1.05;
P = MRF_get_param_dict(P, {'T2<T1'});
look_up       = [P.T1, P.T2];
look_up_names = {'T1', 'T2'};

% T1, T2, T1p
% P.T1.range  = [0.01,  4]; P.T1.factor  = 1.05;
% P.T2.range  = [0.001, 3]; P.T2.factor  = 1.05;
% P.T1p.range = [0.001, 3]; P.T1p.factor = 1.05;
% P = MRF_get_param_dict(P, {'T2<T1', 'T2<T1p', 'T1p<T1'});
% look_up       = [P.T1, P.T2, P.T1p];
% look_up_names = {'T1', 'T2', 'T1p'};

% T1, T2, B1+ correction
% P.T1.range  = [0.01,  4]; P.T1.factor = 1.025;
% P.T2.range  = [0.001, 3]; P.T2.factor = 1.025;
% P.db1.range = [0.8, 1.2]; P.db1.step  = 0.025;
% P = MRF_get_param_dict(P, {'T2<T1'});
% look_up       = [P.T1, P.T2, P.db1];
% look_up_names = {'T1', 'T2', 'db1'};

params_dict.P             = P;
params_dict.look_up       = look_up;
params_dict.look_up_names = look_up_names;
clear P look_up look_up_names;

%% parameters: low-rank reconstruction & matching
params_reco.TestReco           = true;   % do test reco: mixed contrast of all TRs
params_reco.DirectMatching     = false;  % do reco via direct matching
params_reco.DirectMatching_SVD = false;  % do reco via SVD compression of the dictionary before matching
params_reco.LowRank            = true;   % do reco via iterative low-rank reconstruction
params_reco.ROVIR              = true;   % use ROVIR coils for outer FOV artifact suppression
params_reco.rovir_mode         = 'auto'; % 'auto' uses the oversampled region; 'manual' interactive ROI selection
params_reco.rovir_method       = 'auto'; % 'auto' automatic estimation of virtual coils based on thresholding; 'manual' define number of coils
params_reco.rovir_thresh       = 2;      % threshold for automatic coil compression
params_reco.CoilComp           = true;   % additional SVD coil compression after ROVIR
params_reco.NCoils_v           = 8;      % final number of virtual coils after ROVIR and SVD compression
params_reco.ESPIRiT            = true;   % use ESPIRiT or openadapt for calculating cmaps
params_reco.readOS             = 2;      % read oversampling factor
params_reco.NBlocks            = 12;     % number of blocks for pattern matching
params_reco.zero_params.on_off = true;   % parameters for zero interpolation filling
params_reco.zero_params.filter = 'kaiser';
params_reco.zero_params.factor = 2.0;
params_reco.zero_params.radius = 6.0;

% low-rank reconstruction parameters
params_LR = setupParameters_LowrankMRF2D();

%% loop over all sub slices

if params_reco.zero_params.on_off
    Nxy = PULSEQ.FOV.Nxy * params_reco.zero_params.factor;
else
    Nxy = PULSEQ.FOV.Nxy;
end

% init 3D result arrays
M0_Maps = zeros(Nz, Nxy, Nxy);
T1_Maps = zeros(Nz, Nxy, Nxy);
T2_Maps = zeros(Nz, Nxy, Nxy);
IP_Maps = zeros(Nz, Nxy, Nxy);
PC1s    = zeros(Nz, Nxy, Nxy);

for nz = 1:Nz

    % load sub-slice rawdata
    DATA_sub = squeeze(DATA_3D_FFTz(nz, :,:,:));

    % change z axis for isochromat simulation (not relevant for EPG)
    params_dict.z = PULSEQ.FOV.dz / PULSEQ.FOV.Nz * ( linspace(0, 1, params_dict.N_iso)' + nz - PULSEQ.FOV.Nz/2 - 1);    
    
    % simulate dictionary, reconstruction of rawdata & match parameters
    [match, images] = mrf_reco(DATA_sub, NOISE, PULSEQ, params_dict, params_reco, params_LR, ktraj_meas);
    
    % store results in 3D arrays
    M0_Maps(nz,:,:) = match.LR.M0;
    T1_Maps(nz,:,:) = match.LR.T1;
    T2_Maps(nz,:,:) = match.LR.T2;
    IP_Maps(nz,:,:) = match.LR.IP;
    PC1s(nz,:,:)    = squeeze(images.LR(1,:,:));

end
