%% load pulseq backup, read seq file, prepare simulation
clear;
load('Q:/data/Pulseq/Pulseq_Workspace/mgram/260313/260313_0507/backup_260313_0507_workspace.mat')

% find path of the original .seq file
[~, ~, seq_path] = pulseq_get_user_definitions();
seq_path = [seq_path '/Pulseq_Workspace/' PULSEQ.pulseq_user '/' PULSEQ.seq_id(1:6) '/' PULSEQ.seq_id '/' PULSEQ.seq_name '.seq' ];

f0         = PULSEQ.system.B0 * PULSEQ.system.gamma;
[SEQ, SIM] = MRF_read_seq_file( seq_path, ... % .seq file path
                                f0, ...       % larmor frequency
                                [], ...       % adc time stamps
                                [], ...       % soft delays
                                [], ...       % kz part
                                [], ...       % echo mode
                                1e-6, ...     % simulation raster time
                                0);           % plot flag

% load sequence parameters
ktraj = SPI_load_ktraj(PULSEQ);
if isfield(PULSEQ.SPI, 'proj')
    proj_id = PULSEQ.SPI.proj.id;
else
    proj_id = PULSEQ.SPI.phi_id;
end
fov  = PULSEQ.FOV.fov_x;
NR   = PULSEQ.SPI.NR;
Nxy  = PULSEQ.FOV.Nxy;
Nadc = size(ktraj, 3);

%% ---------- design numerical phantom I, simulate Mxy(TR), simulate MRI signals ----------

sim_mode = 'BLOCH'; % 'BLOCH' or 'EPG'
N_iso    = 1000;
s_fac    = 2;

% ---------- step 1: design numerical phantom with physical properties ----------

% define geometry
geo.n_inner = 5;    % number of circular phantom inserts (inner)
geo.n_outer = 9;    % number of circular phantom inserts (outer)
geo.r_phant = 0.9;  % radius of main phantom
geo.r_outer = 0.65; % radius of circle along which outer phantoms lay
geo.r_inner = 0.30; % radius of circle along which inner phantoms lay
geo.r_insrt = 0.1;  % radius of each circular insert
phantom1    = mg_numerical_phantom(Nxy, geo);

% add parameter noise
phantom1.noise_level = 0.01; % [%]

% choose physical properties
[T1, T2]  = NIST_references('3.0T_MnCl2');
temp_ind  = [11; 10; 1; 7; 14; 3; 2; 12; 4; 9; 8; 5; 6; 13];
PD        = ones(1, phantom1.n_phant + 1)';
T1        = [3.0; T1(temp_ind)];
T2        = [1.0; T2(temp_ind)];
T2s       = T2 .* (1+rand(15,1))/2;
phantom1  = mg_phantom_properties(phantom1, PD, T1, T2, T2s);
clear PD temp_ind T1 T2 T2s;

% ---------- step 2: simulate Mxy for all voxels at beginning of all ADCs ----------
params_dict.P = phantom1.p1D;
switch sim_mode
    case 'EPG'
        Mxy = MRF_sim_EPG(MRF_sim_pre(SIM, params_dict.P, [], 'EPG', 0, 0), params_dict.P, 0);
    case 'BLOCH'
        params_dict.z = linspace(-1/2, 1/2, N_iso)' * s_fac * PULSEQ.FOV.dz;
        Mxy = MRF_sim_BLOCH(MRF_sim_pre(SIM, params_dict.P, params_dict.z, 'BLOCH', 1, 0), params_dict.P, params_dict.z, [], 0);
end
clear params_dict;

% ---------- step 3: simulate MRI signals ----------
DATA1 = mg_phantom_sim(Mxy .* repmat(phantom1.p1D.PD,[1,NR]).', phantom1.p2D.mask, phantom1.p2D.cmaps, ktraj, proj_id, phantom1.p2D.T2s, phantom1.p2D.dw0, PULSEQ.SPI.adc.dwell, PULSEQ.FOV.fov_xy);
DATA1 = DATA1 / max(abs(DATA1(:))) * max(abs(Mxy(:)));
clear Mxy;

%% ---------- design numerical phantom II, simulate Mxy(TR), simulate MRI signals ----------

% ---------- step 1: design numerical phantom with physical properties ----------

% use identical geometry
phantom2 = mg_numerical_phantom(Nxy, geo);
clear geo;

% add parameter noise
phantom2.noise_level = 0.01; % [%]

% choose physical properties
[T1, T2]  = NIST_references('3.0T_MnCl2');
temp_ind  = [11; 10; 1; 7; 14; 3; 2; 12; 4; 9; 8; 5; 6; 13];
PD        = [zeros(8,1); linspace(0,1,6)'];
dw0       = [zeros(8,1); f0*3.5*1e-6 * ones(6,1)];
PD        = [0; PD(temp_ind)];
dw0       = [0; dw0(temp_ind)];
T1        = [3.0; T1(temp_ind)] / 3;
T2        = [1.0; T2(temp_ind)] / 3;
T2s       = T2 .* (1+rand(15,1))/2;
phantom2  = mg_phantom_properties(phantom2, PD, T1, T2, T2s, dw0);
clear PD temp_ind T1 T2 T2s;

% ---------- step 2: simulate Mxy for all voxels at beginning of all ADCs ----------
params_dict.P = phantom2.p1D;
switch sim_mode
    case 'EPG'
        Mxy = MRF_sim_EPG(MRF_sim_pre(SIM, params_dict.P, [], 'EPG', 0, 0), params_dict.P, 0);
    case 'BLOCH'
        params_dict.z = linspace(-1/2, 1/2, N_iso)' * s_fac * PULSEQ.FOV.dz;
        Mxy = MRF_sim_BLOCH(MRF_sim_pre(SIM, params_dict.P, params_dict.z, 'BLOCH', 1, 0), params_dict.P, params_dict.z, [], 0);
end
clear params_dict;

% ---------- step 3: simulate MRI signals ----------
DATA2 = mg_phantom_sim(Mxy .* repmat(phantom2.p1D.PD,[1,NR]).', phantom2.p2D.mask, phantom2.p2D.cmaps, ktraj, proj_id, phantom2.p2D.T2s, phantom2.p2D.dw0, PULSEQ.SPI.adc.dwell, PULSEQ.FOV.fov_xy);
DATA2 = DATA2 / max(abs(DATA2(:))) * max(abs(Mxy(:)));
clear Mxy;

%% combine signals and add noise
rng("default")
noise_level_signal = 0.001; % [%]
DATA = DATA1 + DATA2 + (randn(size(DATA1)) + 1i * randn(size(DATA1))) * noise_level_signal / 100;

%% ---------- mrf reco of simulated MRI signals ----------

% parameters: dictionary simulation
params_dict.seq_path    = [];       % full path of .seq file; [] for automatic search in user backups
params_dict.sim_mode    = sim_mode; % 'BLOCH' or 'EPG'
params_dict.comp_energy = 0;        % svd compression energy: 0 for uncompressed dictionary
params_dict.N_iso       = N_iso;    % number of isochromats for bloch simulation
params_dict.s_fac       = s_fac;    % factor for out-of-slice simulation
params_dict.f0          = f0;       % larmor frequency f0: required to simulate rf pulses with ppm freq offsets
params_dict.time_stamps = [];       % adc times stamps: required for correction of trigger delays in cMRF
params_dict.soft_delays = [];       % soft delay user input: required for correction of the acq window in cMRF
params_dict.flag_kz     = [];       % find kz partitions for stacked 3D MRF -> eliminate unnecessary partitions
params_dict.echo_mode   = [];       % echo mode; default: 'spiral_out'

% parameters: dictionary and look-up table
P.T1.range = [0.01,  3.5]; P.T1.factor = 1.05;
P.T2.range = [0.001, 1.5]; P.T2.factor = 1.05;
P = MRF_get_param_dict(P, {'T2<T1'});
look_up       = [P.T1, P.T2];
look_up_names = {'T1', 'T2'};
params_dict.P             = P;
params_dict.look_up       = look_up;
params_dict.look_up_names = look_up_names;
clear P look_up look_up_names;

% parameters: low-rank reconstruction & matching
params_reco.TestReco           = true;   % do test reco: mixed contrast of all TRs
params_reco.DirectMatching     = true;   % do reco via direct matching
params_reco.DirectMatching_SVD = true;   % do reco via SVD compression of the dictionary before matching
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
params_reco.zero_params.on_off = false;  % parameters for zero interpolation filling

% low-rank reconstruction parameters
params_LR = setupParameters_LowrankMRF2D();

%% simulate dictionary, reconstruction of rawdata & match parameters
[match, images] = mrf_reco(DATA, [], PULSEQ, params_dict, params_reco, params_LR, []);

%% deviations: mean +- std [%]

mask_inserts = phantom1.p2D.ind>1;
mask_inserts = imerode(mask_inserts, strel('diamond', 4));

% T1
temp_dev_t1_direct = (match.direct.T1 - phantom1.p2D.T1) ./ phantom1.p2D.T1 * 100;
temp_dev_t1_direct = temp_dev_t1_direct(mask_inserts);
temp_dev_t1_SVD    = (match.SVD.T1 - phantom1.p2D.T1) ./ phantom1.p2D.T1 * 100;
temp_dev_t1_SVD    = temp_dev_t1_SVD(mask_inserts);
temp_dev_t1_LR     = (match.LR.T1 - phantom1.p2D.T1) ./ phantom1.p2D.T1 * 100;
temp_dev_t1_LR     = temp_dev_t1_LR(mask_inserts);
fig1 = figure();
subplot(2,1,1)
hold on
histogram(temp_dev_t1_direct, 'Normalization', 'probability', 'FaceAlpha', 0.5);
histogram(temp_dev_t1_SVD,    'Normalization', 'probability', 'FaceAlpha', 0.5);
histogram(temp_dev_t1_LR,     'Normalization', 'probability', 'FaceAlpha', 0.5);
xlim([-100 100])
ylabel('Relative Frequency')
legend( ['direct: (' num2str(mean(temp_dev_t1_direct), '%.2f') ' ± ' num2str(std(temp_dev_t1_direct), '%.2f') ') %'], ...
        ['SVD: ('    num2str(mean(temp_dev_t1_SVD), '%.2f')    ' ± ' num2str(std(temp_dev_t1_SVD), '%.2f') ') %'], ...
        ['LR: ('     num2str(mean(temp_dev_t1_LR), '%.2f')     ' ± ' num2str(std(temp_dev_t1_LR), '%.2f') ') %'] );
title('T1')
set(gca, 'FontName', 'arial', 'FontSize', 12, 'FontWeight', 'Bold', 'LineWidth', 2)

% T2
temp_dev_t2_direct = (match.direct.T2 - phantom1.p2D.T2) ./ phantom1.p2D.T2 * 100;
temp_dev_t2_direct = temp_dev_t2_direct(mask_inserts);
temp_dev_t2_SVD    = (match.SVD.T2 - phantom1.p2D.T2) ./ phantom1.p2D.T2 * 100;
temp_dev_t2_SVD    = temp_dev_t2_SVD(mask_inserts);
temp_dev_t2_LR     = (match.LR.T2 - phantom1.p2D.T2) ./ phantom1.p2D.T2 * 100;
temp_dev_t2_LR     = temp_dev_t2_LR(mask_inserts);
subplot(2,1,2)
hold on
histogram(temp_dev_t2_direct, 'Normalization', 'probability', 'FaceAlpha', 0.5);
histogram(temp_dev_t2_SVD,    'Normalization', 'probability', 'FaceAlpha', 0.5);
histogram(temp_dev_t2_LR,     'Normalization', 'probability', 'FaceAlpha', 0.5);
xlim([-100 100])
xlabel('Deviation from ground truth [%]')
ylabel('Relative Frequency')
legend( ['direct: (' num2str(mean(temp_dev_t2_direct), '%.2f') ' ± ' num2str(std(temp_dev_t2_direct), '%.2f') ') %'], ...
        ['SVD: ('    num2str(mean(temp_dev_t2_SVD), '%.2f')    ' ± ' num2str(std(temp_dev_t2_SVD), '%.2f') ') %'], ...
        ['LR: ('     num2str(mean(temp_dev_t2_LR), '%.2f')     ' ± ' num2str(std(temp_dev_t2_LR), '%.2f') ') %'] );
title('T2')
set(gca, 'FontName', 'arial', 'FontSize', 12, 'FontWeight', 'Bold', 'LineWidth', 2)

set(gcf, 'WindowState', 'maximized');
% saveas(fig1, [PULSEQ.seq_name '_plot1.svg']);

clear temp_dev_t1_direct temp_dev_t1_SVD temp_dev_t1_LR temp_dev_t2_direct temp_dev_t2_SVD temp_dev_t2_LR;

%% vis match results
t1lims = [0 3000] *1e-3;
t2lims = [0 1100] *1e-3;
t1cmp  = get_cmp('T1', 1000, 1);
t2cmp  = get_cmp('T2', 1000, 1);

fig2 = figure('Name','match results');
ax1 = subplot(3,3,1);
    imagesc(abs(match.direct.M0) .* phantom1.p2D.mask); axis image; axis off; colormap(gca, gray); colorbar;
    title('M0 direct');
ax2 = subplot(3,3,2);
    imagesc(abs(match.SVD.M0) .* phantom1.p2D.mask); axis image; axis off; colormap(gca, gray); colorbar;
    title('M0 SVD');
ax3 = subplot(3,3,3);
    imagesc(abs(match.LR.M0) .* phantom1.p2D.mask); axis image; axis off; colormap(gca, gray); colorbar;
    title('M0 LR');

ax4 = subplot(3,3,4);
    imagesc(match.direct.T1 .* phantom1.p2D.mask, t1lims); axis image; axis off; colormap(gca, t1cmp); colorbar;
    title('T1 direct');
ax5 = subplot(3,3,5);
    imagesc(match.SVD.T1 .* phantom1.p2D.mask, t1lims); axis image; axis off; colormap(gca, t1cmp); colorbar;
    title('T1 SVD');
ax6 = subplot(3,3,6);
    imagesc(match.LR.T1 .* phantom1.p2D.mask, t1lims); axis image; axis off; colormap(gca, t1cmp); colorbar;
    title('T1 LR');

ax7 = subplot(3,3,7);
    imagesc(match.direct.T2 .* phantom1.p2D.mask, t2lims); axis image; axis off; colormap(gca, t2cmp); colorbar;
    title('T2 direct');
ax8 = subplot(3,3,8);
    imagesc(match.SVD.T2 .* phantom1.p2D.mask, t2lims); axis image; axis off; colormap(gca, t2cmp); colorbar;
    title('T2 SVD');
ax9 = subplot(3,3,9);
    imagesc(match.LR.T2 .* phantom1.p2D.mask, t2lims); axis image; axis off; colormap(gca, t2cmp); colorbar;
    title('T2 LR');

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9]);
clear ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9;

set(gcf, 'WindowState', 'maximized');
% saveas(fig2, [PULSEQ.seq_name '_plot2.svg']);

%% vis difference vs ground truth

dlims = [-20 20]; % [%]
dcmp  = [[linspace(0,1,1000), linspace(1,1,1000)]; ...
        [ linspace(0,1,1000), linspace(1,0,1000)]; ...
        [ linspace(1,1,1000), linspace(1,0,1000)]]';

dT1_direct = (match.direct.T1 - phantom1.p2D.T1) ./ phantom1.p2D.T1 *100;
dT1_SVD    = (match.SVD.T1    - phantom1.p2D.T1) ./ phantom1.p2D.T1 *100;
dT1_LR     = (match.LR.T1     - phantom1.p2D.T1) ./ phantom1.p2D.T1 *100;
dT2_direct = (match.direct.T2 - phantom1.p2D.T2) ./ phantom1.p2D.T2 *100;
dT2_SVD    = (match.SVD.T2    - phantom1.p2D.T2) ./ phantom1.p2D.T2 *100;
dT2_LR     = (match.LR.T2     - phantom1.p2D.T2) ./ phantom1.p2D.T2 *100;

dT1_direct(phantom1.p2D.mask==0) = 0;
dT1_SVD(phantom1.p2D.mask==0)    = 0;
dT1_LR(phantom1.p2D.mask==0)     = 0;
dT2_direct(phantom1.p2D.mask==0) = 0;
dT2_SVD(phantom1.p2D.mask==0)    = 0;
dT2_LR(phantom1.p2D.mask==0)     = 0;

fig3 = figure('Name','match deviations');
ax1 = subplot(2,3,1);
    imagesc(dT1_direct, dlims); axis image; axis off; colormap(dcmp); colorbar;
    title('T1 direct');

ax2 = subplot(2,3,2);
    imagesc(dT1_SVD, dlims); axis image; axis off; colormap(dcmp); colorbar;
    title('T1 SVD');

ax3 =  subplot(2,3,3);
    imagesc(dT1_LR, dlims); axis image; axis off; colormap(dcmp); colorbar;
    title('T1 LR');

ax4 = subplot(2,3,4);
    imagesc(dT2_direct, dlims); axis image; axis off; colormap(dcmp); colorbar;
    title('T2 direct');

ax5 = subplot(2,3,5);
    imagesc(dT2_SVD, dlims); axis image; axis off; colormap(dcmp); colorbar;
    title('T2 SVD');

ax6 = subplot(2,3,6);
    imagesc(dT2_LR, dlims); axis image; axis off; colormap(dcmp); colorbar;
    title('T2 LR'); 

linkaxes([ax1 ax2 ax3 ax4 ax4 ax5 ax6]);
clear ax1 ax2 ax3 ax4 ax5 ax6;

set(gcf, 'WindowState', 'maximized');
% saveas(fig3, [PULSEQ.seq_name '_plot3.svg']);

%% vis difference vs ground truth

for j=1:phantom1.n_phant
    temp = match.SVD.T1((phantom1.p2D.ind==j+1).*mask_inserts==1);
    t1_vals_SVD(j,1) = mean(temp);
    t1_vals_SVD(j,2) = std(temp);
    temp = match.SVD.T2((phantom1.p2D.ind==j+1).*mask_inserts==1);
    t2_vals_SVD(j,1) = mean(temp);
    t2_vals_SVD(j,2) = std(temp);
end
for j=1:phantom1.n_phant
    temp = match.LR.T1((phantom1.p2D.ind==j+1).*mask_inserts==1);
    t1_vals_LR(j,1) = mean(temp);
    t1_vals_LR(j,2) = std(temp);
    temp = match.LR.T2((phantom1.p2D.ind==j+1).*mask_inserts==1);
    t2_vals_LR(j,1) = mean(temp);
    t2_vals_LR(j,2) = std(temp);
end
clear temp;

fig4 = figure();
ax1 = subplot(2,2,1);
imagesc(dT1_LR, dlims); axis image; axis off; colormap(dcmp); colorbar; title('T1 deviation [%]');

ax2 = subplot(2,2,2);
imagesc(dT2_LR, dlims); axis image; axis off; colormap(dcmp); colorbar; title('T2 deviation [%]');

ax3 = subplot(2,2,3);
hold on;
plot([1e-4, 10], [1e-4, 10], 'k-', 'LineWidth', 2);
plot(phantom1.p.T1(2:phantom1.n_phant+1), t1_vals_SVD(:,1), 'r.', 'MarkerSize', 25);
plot(phantom1.p.T1(2:phantom1.n_phant+1), t1_vals_LR(:,1),  'b.', 'MarkerSize', 25);
xlim([0.06, 4]); ylim([0.06, 4]); xlabel('T1 ground truth [s]'); ylabel('T1 mrf [s]'); ax3.XScale='log'; ax3.YScale='log'; grid on; axis square;
set(gca, 'Fontname', 'arial', 'Fontweight', 'bold', 'Fontsize', 12, 'LineWidth', 2);

ax4 = subplot(2,2,4);
hold on;
plot([1e-4, 10], [1e-4, 10], 'k-', 'LineWidth', 2);
plot(phantom1.p.T2(2:phantom1.n_phant+1), t2_vals_SVD(:,1), 'r.', 'MarkerSize', 25);
plot(phantom1.p.T2(2:phantom1.n_phant+1), t2_vals_LR(:,1),  'b.', 'MarkerSize', 25);
plot([1e-4, 10], [1e-4, 10], 'k--', 'LineWidth', 2);
xlim([0.005, 1.5]); ylim([0.005, 1.5]); xlabel('T2 ground truth [s]'); ylabel('T2 mrf [s]'); ax4.XScale='log'; ax4.YScale='log'; grid on; axis square;
set(gca, 'Fontname', 'arial', 'Fontweight', 'bold', 'Fontsize', 12, 'LineWidth', 2);
clear ax1 ax2 ax3 ax4;

set(gcf, 'WindowState', 'maximized');
% saveas(fig4, [PULSEQ.seq_name '_plot4.svg']);