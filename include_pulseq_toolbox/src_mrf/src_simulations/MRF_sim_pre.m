function SIM = MRF_sim_pre(SIM, P, z, sim_mode, flag_plot_slice, flag_plot_adiabatic)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% this function compresses the operator lists for fast bloch or epg simulations.
% global pulses, slice selective pulses and adibatic pulses are compressed
% to instantaneous events. for correction, the slice profiles and
% relaxation loss during adibatic pulses are pre-simulated and stored in
% the SPROF array or the T1/T2 dependent transition matrix TM
% the final outputs can be simulated with MRF_sim_BLOCH() or MRF_sim_EPG()

% ----- input: -----
% SIM -> generated with MRF_read_seq_file()
%   .ID:  compressed ID array for selecting the correct Bloch operator
%   .RF:  compressed complex RF waveforms [rad/s]
%   .GZ:  compressed z gradient moments [rad/m] or waveforms [rad/s/m]
%   .DB:  compressed b-value of the diffusion encoding tensor [s/mm^2]
%   .DT:  compressed variable time steps [s]
%   .PHI: adc reciever phases [rad]
% P contains the T1 and T2 values for the pre-simulations
%   .T1: N_dict x 1
%   .T2: N_dict x 1
% z: N_iso  x 1; z-axis for pre-simulation of slice profiles
% sim_mode: 'BLOCH' or 'EPG'

% ----- output: -----
% SIM 
%   .ID, .RF, .GZ, .DB, .DT, .PHI -> compressed, for faster bloch or epg simulations
%   .SPROF: N_iso x 270 -> slice profile array for flip anggles 1...270°
%   .TM:    N_dict x N_pulses x 4x4 (or 3x3 for EPG) -> transition matrix for all adiabaitc pulses and T1/T2

% ----- step I -----
% compress ID, RF, GZ, DB, DT
% RF waveforms will be replaced by flip angles, phases or adiabatic pulse IDs
% ID==2:  global standard  -> calc flip angles and phases
% ID==3:  global adiabatic -> extract all waveforms and generate IDs
% ID==4:  global spin-lock -> already compressed
% ID==5:  slice selective excitation -> extract averaged normalized waveform, calc flip angles and phases
% ID==6:  slice selective refocusing -> extract averaged normalized waveform, calc flip angles and phases
% ID==12: adiabatic spin-lock -> save duration for simulating adiabatic T1p

% ----- step II -----
% pre-simulation for RF pulses
% wave_adia, wave_exc, wave_rfc will be used for pre-simulations
% wave_adia -> T1/T2 depending calculation of the transition matrix
% wave_exc  -> slice profile simulation for 1° 2° ... 135°
% wave_rfc  -> slice profile simulation for 136° 137° ... 270°

% defaults
if nargin<6
    flag_plot_adiabatic = [];
end
if nargin<5
    flag_plot_slice = [];
end
if nargin<4
    sim_mode = [];
end
if isempty(flag_plot_adiabatic)
    flag_plot_adiabatic = 0;
end
if isempty(flag_plot_slice)
    flag_plot_slice = 0;
end
if isempty(sim_mode)
    sim_mode = 'BLOCH';
end

ID  = SIM.ID;
RF  = SIM.RF;
GZ  = SIM.GZ;
DB  = SIM.DB;
DT  = SIM.DT;
PHI = SIM.PHI;
seq_name = SIM.seq_name;
clear SIM;

%% -------------------------- step I: compress --------------------------

% grad raster time
dt = unique(DT);
dt(dt==0) = [];
dt = dt(1);

%% search for global standard pulses: calc FAs and phases
temp_id    = ID==2;
temp_id    = [0; temp_id; 0];
temp_id    = diff(temp_id);
temp_start = find(temp_id==1);
temp_end   = find(temp_id==-1)-1;
for j=1:numel(temp_start)
    temp_wave   = RF(temp_start(j):temp_end(j));
    temp_fa     = abs(sum(temp_wave)) * dt;
    temp_center = find(abs(temp_wave)==max(abs(temp_wave)));
    temp_center = floor(mean(temp_center));    
    temp_phase  = angle(temp_wave(temp_center));
    temp_fa     = temp_fa * exp(1i*temp_phase);    
    temp_dur1   = temp_center * dt;
    temp_dur2   = (numel(temp_wave)-temp_center) * dt;    
    RF(temp_start(j)+0) = 0;
    RF(temp_start(j)+1) = temp_fa; % store flip angle and phase
    RF(temp_start(j)+2) = 0;    
    DT(temp_start(j)+0) = temp_dur1; % store duration for free relaxation & dephasing
    DT(temp_start(j)+1) = 0;
    DT(temp_start(j)+2) = temp_dur2; % store duration for free relaxation & dephasing        
    ID(temp_start(j)+0) = 1;  % free relaxation & dephasing
    ID(temp_start(j)+1) = 2;  % instantaneous RF excitation
    ID(temp_start(j)+2) = 1;  % free relaxation & dephasing    
    ID(temp_start(j)+3:temp_end(j)) = -1; % remove IDs, compress simulation
    if numel(temp_wave)<4
        error('rf pulse is too short for compression');
    end
end
clear j temp_id temp_start temp_end temp_n temp_wave temp_fa temp_center temp_phase temp_dur1 temp_dur2;

%% search for adiabatic pulses: extract waveforms and generate IDs
temp_id    = ID==3;
temp_id    = [0; temp_id; 0];
temp_id    = diff(temp_id);
temp_start = find(temp_id==1);
temp_end   = find(temp_id==-1)-1;
wave_adia  = struct();
count      = 1;
for j=1:numel(temp_start)
    temp_wave = RF(temp_start(j):temp_end(j));
    temp_hash = pulseq_get_wave_hash([real(temp_wave); imag(temp_wave)]);
    if ~isfield(wave_adia, 'hash')
        temp_no = 1;
        wave_adia(count,1).wave = temp_wave;
        wave_adia(count,1).hash = temp_hash;
        count = count + 1;
    else
        temp_check = 0;
        for k=1:numel(wave_adia)
            if strcmp(temp_hash, wave_adia(k).hash)
                temp_check = 1;
                temp_no = k;
            end
        end
        if temp_check==0
            temp_no = count;
            wave_adia(count,1).wave = temp_wave;
            wave_adia(count,1).hash = temp_hash;
            count = count + 1;
        end
    end
    RF(temp_start(j)) = temp_no;
    DT(temp_start(j)) = sum(DT(temp_start(j):temp_end(j)));
    ID(temp_start(j)+1:temp_end(j)) = -1;    
end
clear j k temp_id temp_start temp_end temp_n temp_wave temp_hash temp_check temp_no count;

%% search for slice selective excitation pulses: extract waveforms and calc FAs and phases
temp_id    = ID==5;
temp_id    = [0; temp_id; 0];
temp_id    = diff(temp_id);
temp_start = find(temp_id==1);
temp_end   = find(temp_id==-1)-1;
temp_n     = numel(temp_start);
wave_exc   = 0;
exc_exist  = 0;
for j=1:numel(temp_start)
    temp_wave   = RF(temp_start(j):temp_end(j));
    temp_fa     = abs(sum(temp_wave)) * dt;
    temp_center = find(abs(temp_wave)==max(abs(temp_wave)));
    temp_center = floor(mean(temp_center));    
    temp_phase  = angle(temp_wave(temp_center));
    temp_fa     = temp_fa * exp(1i*temp_phase);    
    temp_wave   = temp_wave * exp(-1i*temp_phase);
    wave_exc    = wave_exc + temp_wave / max(abs(temp_wave)) / temp_n;
    temp_dur1   = temp_center * dt;
    temp_dur2   = (numel(temp_wave)-temp_center) * dt;
    temp_gz     = mean(GZ(temp_start(j):temp_end(j)));
    RF(temp_start(j)+0) = 0;
    RF(temp_start(j)+1) = 0;
    RF(temp_start(j)+2) = temp_fa; % store flip angle and phase
    RF(temp_start(j)+3) = 0;
    RF(temp_start(j)+4) = 0;
    DT(temp_start(j)+0) = temp_dur1; % store duration for free relaxation & dephasing
    DT(temp_start(j)+1) = 0;
    DT(temp_start(j)+2) = 0;
    DT(temp_start(j)+3) = 0;
    DT(temp_start(j)+4) = temp_dur2; % store duration for free relaxation & dephasing
    GZ(temp_start(j)+0) = 0;
    GZ(temp_start(j)+1) = temp_gz * temp_dur1; % store gradient moment
    GZ(temp_start(j)+2) = temp_gz;             % store gradient strength
    GZ(temp_start(j)+3) = temp_gz * temp_dur1; % store gradient moment
    GZ(temp_start(j)+4) = 0;
    ID(temp_start(j)+0) = 1;  % free relaxation & dephasing
    ID(temp_start(j)+1) = 10; % z gradient dephasing
    ID(temp_start(j)+2) = 5;  % instantaneous RF excitation
    ID(temp_start(j)+3) = 10; % z gradient dephasing
    ID(temp_start(j)+4) = 1;  % free relaxation & dephasing
    ID(temp_start(j)+5:temp_end(j)) = -1; % remove IDs, compress simulation
    if temp_fa*180/pi > 135
        warning('excitation pulse > 135°');
    end
    if numel(temp_wave)<6
        error('rf pulse is too short for compression');
    end
    exc_exist = 1;
end
clear j temp_id temp_start temp_end temp_n temp_wave temp_fa temp_center temp_phase temp_dur1 temp_dur2 temp_gz;

%% search for slice selective refocusing pulses: extract waveforms and calc FAs and phases
temp_id    = ID==6;
temp_id    = [0; temp_id; 0];
temp_id    = diff(temp_id);
temp_start = find(temp_id==1);
temp_end   = find(temp_id==-1)-1;
temp_n     = numel(temp_start);
wave_rfc   = 0;
rfc_exist  = 0;
for j=1:numel(temp_start)
    temp_wave   = RF(temp_start(j):temp_end(j));
    temp_fa     = abs(sum(temp_wave)) * dt;
    temp_center = find(abs(temp_wave)==max(abs(temp_wave)));
    temp_center = floor(mean(temp_center));    
    temp_phase  = angle(temp_wave(temp_center));
    temp_fa     = temp_fa * exp(1i*temp_phase);    
    temp_wave   = temp_wave * exp(-1i*temp_phase);
    wave_rfc    = wave_rfc + temp_wave / max(abs(temp_wave)) / temp_n;
    temp_dur1   = temp_center * dt;
    temp_dur2   = (numel(temp_wave)-temp_center) * dt;
    temp_gz     = mean(GZ(temp_start(j):temp_end(j)));
    RF(temp_start(j)+0) = 0;
    RF(temp_start(j)+1) = 0;
    RF(temp_start(j)+2) = temp_fa; % store flip angle and phase
    RF(temp_start(j)+3) = 0;
    RF(temp_start(j)+4) = 0;
    DT(temp_start(j)+0) = temp_dur1; % store duration for free relaxation & dephasing
    DT(temp_start(j)+1) = 0;
    DT(temp_start(j)+2) = 0;
    DT(temp_start(j)+3) = 0;
    DT(temp_start(j)+4) = temp_dur2; % store duration for free relaxation & dephasing
    GZ(temp_start(j)+0) = 0;
    GZ(temp_start(j)+1) = temp_gz * temp_dur1; % store gradient moment
    GZ(temp_start(j)+2) = temp_gz;             % store gradient strength
    GZ(temp_start(j)+3) = temp_gz * temp_dur1; % store gradient moment
    GZ(temp_start(j)+4) = 0;
    ID(temp_start(j)+0) = 1;  % free relaxation & dephasing
    ID(temp_start(j)+1) = 10; % z gradient dephasing
    ID(temp_start(j)+2) = 5;  % instantaneous RF refocusing
    ID(temp_start(j)+3) = 10; % z gradient dephasing
    ID(temp_start(j)+4) = 1;  % free relaxation & dephasing
    ID(temp_start(j)+5:temp_end(j)) = -1; % remove IDs, compress simulation
    if temp_fa*180/pi < 135
        warning('excitation pulse > 135°');
    end
    if numel(temp_wave)<6
        error('rf pulse is too short for compression');
    end    
    rfc_exist = 1;
end
clear j temp_id temp_start temp_end temp_n temp_wave temp_fa temp_center temp_phase;

%% search for adiabatic spin-lock pulses: save duration
temp_id    = ID==12;
temp_id    = [0; temp_id; 0];
temp_id    = diff(temp_id);
temp_start = find(temp_id==1);
temp_end   = find(temp_id==-1)-1;
for j=1:numel(temp_start)
    temp_wave         = RF(temp_start(j):temp_end(j));    
    temp_dur          = numel(temp_wave) * dt;
    RF(temp_start(j)) = 0;
    DT(temp_start(j)) = temp_dur; % store duration for simulating adiabatic T1p
    ID(temp_start(j)) = 12;  % adiabatic spin-lock
    ID(temp_start(j)+1:temp_end(j)) = -1; % remove IDs, compress simulation
end
clear j temp_id temp_start temp_end temp_wave temp_dur;

%% compress
RF(ID==-1) = [];
GZ(ID==-1) = [];
DB(ID==-1) = [];
DT(ID==-1) = [];
ID(ID==-1) = [];

%% join relaxation, spoiler, crusher and diffusion weighting

% init joined arrays
ID_join    = zeros(size(ID)) - 1;
RF_join    = zeros(size(RF));
GZ_join    = zeros(size(GZ));
DB_join    = zeros(size(DB));
DT_join    = zeros(size(DT));
ID_join(1) = ID(1);
RF_join(1) = RF(1);
GZ_join(1) = GZ(1);
DB_join(1) = DB(1);
DT_join(1) = DT(1);
j_start    = 1;
idx        = 2;
temp_rf_adc_gx_gy = ismember(ID, [0 2 3 4 5 6 7 8 9 12]);

while j_start<numel(ID)
    j_end = j_start + 1;
    while temp_rf_adc_gx_gy(j_end)==0
        j_end = j_end + 1;
        if j_end>numel(ID)
            j_end = numel(ID);
            break;
        end
    end
    if j_end > j_start+1
        temp_ind = j_start+1 : j_end-1;
        if abs(sum(DT(temp_ind))) > 0
            ID_join(idx) = 1;
            RF_join(idx) = 0;
            GZ_join(idx) = 0;
            DB_join(idx) = 0;
            DT_join(idx) = sum(DT(temp_ind));
            idx = idx + 1;
        end
        if abs(sum(GZ(temp_ind))) > 0
            ID_join(idx) = 10;
            RF_join(idx) = 0;
            GZ_join(idx) = sum(GZ(temp_ind));
            DB_join(idx) = 0;
            DT_join(idx) = 0;
            idx = idx + 1;
        end        
        if abs(sum(DB(temp_ind))) > 0
            ID_join(idx) = 11;
            RF_join(idx) = 0;
            GZ_join(idx) = 0;
            DB_join(idx) = sum(DB(temp_ind));
            DT_join(idx) = 0;
            idx = idx + 1;
        end
    end
    ID_join(idx) = ID(j_end);
    RF_join(idx) = RF(j_end);
    GZ_join(idx) = GZ(j_end);
    DB_join(idx) = DB(j_end);
    DT_join(idx) = DT(j_end);
    idx = idx + 1;
    j_start = j_end ;
end

% remove unnecessary entries
RF_join(ID_join==-1) = [];
GZ_join(ID_join==-1) = [];
DB_join(ID_join==-1) = [];
DT_join(ID_join==-1) = [];
ID_join(ID_join==-1) = [];

ID = int32(ID_join);
RF = RF_join;
GZ = GZ_join;
DB = DB_join;
DT = DT_join;

clear ID_join RF_join GZ_join DB_join DT_join j_start j_end temp_ind temp_rf_adc_gx_gy;

%% eliminate nearly zero z gradients; abs(moment) < 1 (this means: 1/m = 360°/m = 3.6°/cm) 
temp_del = (abs(GZ)<1) .* (ID==10);
ID(temp_del==1) = [];
RF(temp_del==1) = [];
GZ(temp_del==1) = [];
DB(temp_del==1) = [];
DT(temp_del==1) = [];
clear temp_del;

%% EPG mode: convert to unit gradients
if strcmp(sim_mode, 'EPG')
    GZ(ID~=10) = 0;
    gz_vals    = GZ(ID==10);
    gz_unit    = mode(gz_vals);
    gz_unit    = gz_unit(1);
    GZ         = round(GZ/gz_unit);
end

%% -------------------------- step II: pre-sims --------------------------

T1 = P.T1;
T2 = P.T2;

% efficiency and relaxation loss of adiabatic pulses
if isfield(wave_adia, 'wave')
    TM = pre_sim_adiabatic_pulses(wave_adia, T1, T2, dt, sim_mode);
else
    TM = [];
end

% slice profile correction: excitation and refocusing pulses
if strcmp(sim_mode, 'BLOCH')
    if exc_exist==0 && rfc_exist==0
        SPROF = [];
    elseif exc_exist==1 && rfc_exist==0
        SPROF = presim_slice_profile_exc(ID, GZ, wave_exc, z, dt);
    elseif exc_exist==0 && rfc_exist==1
        SPROF = presim_slice_profile_rfc(ID, GZ, wave_exc, z, dt);
        SPROF = [zeros(N_iso, 135), SPROF];
    else
        SPROF_exc = presim_slice_profile_exc(ID, GZ, wave_exc, z, dt);
        SPROF_rfc = presim_slice_profile_rfc(ID, GZ, wave_exc, z, dt);
        SPROF     = [SPROF_exc, SPROF_rfc];
    end
else
    SPROF = [];
    flag_plot_slice = 0;
end

%% vis: slice profiles
if flag_plot_slice==1
    figure()
    subplot(1,2,1)
    plot((1:numel(wave_exc))'*dt*1e6, abs(wave_exc), 'k-', 'LineWidth', 2)
    xlabel('time t [us]')
    ylabel('norm. pulse amplitude')
    title('rf pulse waveform')
    set(gca, 'FontName', 'arial', 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 2)

    subplot(1,2,2)
    hold on
    for j=1:size(SPROF,2)
        plot(z *1e3, real(SPROF(:,j)), 'b-', 'LineWidth', 1)
        plot(z *1e3, imag(SPROF(:,j)), 'r-', 'LineWidth', 1)
    end
    yline( 1, 'k--')
    yline( 0, 'k--')
    legend('real()', 'imag()')
    xlabel('slice position [mm]')
    ylabel('normalized profile')
    title('simulated slice profiles 1...135°')
    set(gca, 'FontName', 'arial', 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 2)
    hold off
end

%% vis: transition matrix of adiabatic pulses
if flag_plot_adiabatic==1
if isfield(wave_adia, 'wave')    

    for j=1:numel(wave_adia)
        if strcmp(sim_mode, 'BLOCH')
            figure()
            subplot(3,7,[1 2 3 8 9 10 15 16 17])
            hold on
            plot((1:numel(wave_adia(j).wave))*dt/1e-6, real(wave_adia(j).wave/2/pi), '-b', 'LineWidth', 3)
            plot((1:numel(wave_adia(j).wave))*dt/1e-6, imag(wave_adia(j).wave/2/pi), '-r', 'LineWidth', 3)
            xlabel('pulse duration [us]')
            ylabel('real/imag RF [Hz]')
            set(gca, 'FontName', 'Arial', 'FontWeight','Bold')
            subplot(3,7,4);  plot(T2, squeeze(TM(:,j,1,1)), 'k.'); axis square; ylim([-1 1]);
            subplot(3,7,5);  plot(T2, squeeze(TM(:,j,1,2)), 'k.'); axis square; ylim([-1 1]);
            subplot(3,7,6);  plot(T2, squeeze(TM(:,j,1,3)), 'k.'); axis square; ylim([-1 1]);
            subplot(3,7,7);  plot(T2, squeeze(TM(:,j,1,4)), 'k.'); axis square; ylim([-1 1]);
            subplot(3,7,11); plot(T2, squeeze(TM(:,j,2,1)), 'k.'); axis square; ylim([-1 1]);
            subplot(3,7,12); plot(T2, squeeze(TM(:,j,2,2)), 'k.'); axis square; ylim([-1 1]);
            subplot(3,7,13); plot(T2, squeeze(TM(:,j,2,3)), 'k.'); axis square; ylim([-1 1]);
            subplot(3,7,14); plot(T2, squeeze(TM(:,j,2,4)), 'k.'); axis square; ylim([-1 1]);
            subplot(3,7,18); plot(T2, squeeze(TM(:,j,3,1)), 'k.'); axis square; ylim([-1 1]);
            subplot(3,7,19); plot(T2, squeeze(TM(:,j,3,2)), 'k.'); axis square; ylim([-1 1]);
            subplot(3,7,20); plot(T2, squeeze(TM(:,j,3,3)), 'k.'); axis square; ylim([-1 1]);
            subplot(3,7,21); plot(T2, squeeze(TM(:,j,3,4)), 'k.'); axis square; ylim([-1 1]);
            sgtitle(['BLOCH: ' wave_adia(j).hash '   x-axis: T2 [s]'], 'interpreter', 'none')
        end
        
        if strcmp(sim_mode, 'EPG')
            figure()
            subplot(3,6,[1 2 3 7 8 9 13 14 15])
            hold on
            plot((1:numel(wave_adia(j).wave))*dt/1e-6, real(wave_adia(j).wave/2/pi), '-b', 'LineWidth', 3)
            plot((1:numel(wave_adia(j).wave))*dt/1e-6, imag(wave_adia(j).wave/2/pi), '-r', 'LineWidth', 3)
            xlabel('pulse duration [us]')
            ylabel('real/imag RF [Hz]')
            set(gca, 'FontName', 'Arial', 'FontWeight','Bold')
            subplot(3,6,4);  hold on; plot(T2, squeeze(real(TM(:,j,1,1))), 'b.'); plot(T2, squeeze(imag(TM(:,j,1,1))), 'r.'); axis square; ylim([-1 1]);
            subplot(3,6,5);  hold on; plot(T2, squeeze(real(TM(:,j,1,2))), 'b.'); plot(T2, squeeze(imag(TM(:,j,1,2))), 'r.'); axis square; ylim([-1 1]);
            subplot(3,6,6);  hold on; plot(T2, squeeze(real(TM(:,j,1,3))), 'b.'); plot(T2, squeeze(imag(TM(:,j,1,3))), 'r.'); axis square; ylim([-1 1]);
            subplot(3,6,10); hold on; plot(T2, squeeze(real(TM(:,j,2,1))), 'b.'); plot(T2, squeeze(imag(TM(:,j,2,1))), 'r.'); axis square; ylim([-1 1]);
            subplot(3,6,11); hold on; plot(T2, squeeze(real(TM(:,j,2,2))), 'b.'); plot(T2, squeeze(imag(TM(:,j,2,2))), 'r.'); axis square; ylim([-1 1]);
            subplot(3,6,12); hold on; plot(T2, squeeze(real(TM(:,j,2,3))), 'b.'); plot(T2, squeeze(imag(TM(:,j,2,3))), 'r.'); axis square; ylim([-1 1]);
            subplot(3,6,16); hold on; plot(T2, squeeze(real(TM(:,j,3,1))), 'b.'); plot(T2, squeeze(imag(TM(:,j,3,1))), 'r.'); axis square; ylim([-1 1]);
            subplot(3,6,17); hold on; plot(T2, squeeze(real(TM(:,j,3,2))), 'b.'); plot(T2, squeeze(imag(TM(:,j,3,2))), 'r.'); axis square; ylim([-1 1]);
            subplot(3,6,18); hold on; plot(T2, squeeze(real(TM(:,j,3,3))), 'b.'); plot(T2, squeeze(imag(TM(:,j,3,3))), 'r.'); axis square; ylim([-1 1]);    
            sgtitle(['EPG: ' wave_adia(j).hash '   x-axis: T2 [s]'], 'interpreter', 'none')
        end
    end

end
end

%% output SIM struct
SIM.ID  = ID;
SIM.RF  = RF;
SIM.GZ  = GZ;
SIM.DB  = DB;
SIM.DT  = DT;
SIM.PHI = PHI;
SIM.TM  = TM;
if strcmp(sim_mode, 'BLOCH')
    SIM.SPROF = SPROF;
end
SIM.seq_name = seq_name;

end

%% --------------------------------------------------------------------
function TM = pre_sim_adiabatic_pulses(wave_adia, T1, T2, dt, sim_mode)

    % ----- input -----
    % wave_adia: contains all adiabatic rf objects of the sequence
    % T1:        [s] T1 relaxation times
    % T2:        [s] T2 relaxation times
    % dt:        [s] raster time
    % sim_mode:  'BLOCH' or 'EPG'

    % ----- output -----
    % TM: transition matrix for different pulses and relaxation times
    % case BLOCH: N_dict x N_pulses x 4x4
    % case EPG:   N_dict x N_pulses x 3x3    
   
    N_pulses = numel(wave_adia);
    N_dict   = numel(T1);
    if strcmp(sim_mode, 'BLOCH')
        TM = zeros(N_dict, N_pulses, 4, 4);
    elseif strcmp(sim_mode, 'EPG')
        TM = zeros(N_dict, N_pulses, 3, 3);
    else
        error('sim_mode: choose BLOCH or EPG');
    end
       
    % define T1/T2 look-up table
    T1_lookup  = 0.01  * 1.01.^(0:1000); % 1% steps
    T2_lookup  = 0.001 * 1.01.^(0:1000); % 1% steps
    T1_lookup(T1_lookup>6) = [];
    T2_lookup(T2_lookup>4) = [];
    [T1_lookup, T2_lookup] = ndgrid(T1_lookup, T2_lookup);
    T1_lookup  = T1_lookup(:);
    T2_lookup  = T2_lookup(:);
    tempdel = (T2_lookup./T1_lookup) > 1;
    T1_lookup(tempdel) = [];
    T2_lookup(tempdel) = [];
    n_lookup = numel(T1_lookup);
    clear tempdel
    
    % find nearest neighbour for the actual dictionary
    if exist('knnsearch.m', 'file')==2
        ind_lookup = knnsearch([T1_lookup(:), T2_lookup(:)], [T1(:), T2(:)]); % Statistics Toolbox
    else
        ind_lookup = zeros(length(T1), 1);
        parfor j = 1:length(T1)
            [~, ind_lookup(j)] = min(sqrt((T1_lookup - T1(j)).^2 + (T2_lookup - T2(j)).^2)); % slow
        end
    end
       
    for j = 1:N_pulses
    
        % load transition matrix from .mat file
        temp_name = ['transition_matrix_' num2str(round(numel(wave_adia(j).wave)*dt/1e-6)) 'us_' num2str(round(max(abs(wave_adia(j).wave))/2/pi)) 'Hz_' wave_adia(j).hash];
        temp_path = [pulseq_get_path('MRF_find_adia_eff') temp_name '.mat'];    
        if isfile(temp_path)
             temp = load(temp_path);
             if all(T1_lookup==temp.T1_lookup) && all(T2_lookup==temp.T2_lookup)
                 TM_BLOCH = temp.TM_BLOCH;
                 TM_EPG   = temp.TM_EPG;
             end
             clear temp;
        end
    
        % start new simulation and save to .mat file
        if ~(exist('TM_BLOCH', 'var') && exist('TM_EPG', 'var'))
            TM_BLOCH = zeros(n_lookup, 4, 4);
            TM_EPG   = zeros(n_lookup, 3, 3);
    
            % calculate transition matrix for look-up table
            temp_w1x = real(wave_adia(j).wave);
            temp_w1y = imag(wave_adia(j).wave);
            temp_n   = numel(temp_w1x);
            parfor k=1:n_lookup
                [~, temp_trans] = mex_BLOCH_expm(zeros(4,1), temp_w1x, temp_w1y, zeros(temp_n,1), 1/T1_lookup(k), 1/T2_lookup(k), dt);
                temp_trans      = temp_trans .* [1, 1, 1, 1; 1, 1, 1, 1; 1, 1, 1, 1; 0, 0, 0, 1]; % M0=1, d/dt M0=0
                TM_BLOCH(k,:,:) = temp_trans;
            end
            parfor k=1:n_lookup
                temp_trans    = mex_EPG_rf_relax(temp_w1x, temp_w1y, 1/T1_lookup(k), 1/T2_lookup(k), dt);
                TM_EPG(k,:,:) = temp_trans;
            end
    
            % save in mat file
            waveform = wave_adia(j).wave;
            save(temp_path, 'TM_BLOCH', 'TM_EPG', 'T1_lookup', 'T2_lookup', 'waveform', 'dt');            
        end        
   
        % interpolate via lookup table
        if strcmp(sim_mode, 'BLOCH')
            TM(:,j,:,:) = TM_BLOCH(ind_lookup,:,:);
        else
            TM(:,j,:,:) = TM_EPG(ind_lookup,:,:);            
        end

        clear temp_name temp_path TM_BLOCH TM_EPG temp_w1x temp_w1y temp_n waveform;

    end    

end

%% -------------------------------------------------------------------
function SPROF_exc = presim_slice_profile_exc(ID, GZ, wave_exc, z, dt)

    % check equal gradient strength
    if numel(unique(GZ(ID==5)))>1
        error('different gradient strength detected!');
    else
        gz = unique(GZ(ID==5));
    end
    N_iso = numel(z);
    
    % search for slice profile simulation backup
    SProf_hash = ['slice_prof_exc_' ...
                  num2str(numel(wave_exc)) 'us_' ...
                  num2str(gz/1e3) 'kHzm_' ...
                  num2str(N_iso) 'Isos_' ...
                  num2str(min(z)*1e3) '_' num2str(max(z)*1e3) 'mm_' ...
                  pulseq_get_wave_hash(round(abs(wave_exc),3))];
        
    SProf_hash  = [strrep(SProf_hash, '.', ',') '.mat'];
    SProf_path  = pulseq_get_path('MRF_find_slice_profile');
    SProf_path  = [SProf_path SProf_hash];
    SProf_exist = isfile(SProf_path);
    
    if SProf_exist==1
        load(SProf_path);
    else
        % calculate rf waveforms: 1° ... 135°
        FAs = (1:1:135)' *pi/180;
        
        % bloch simulation for isochromats        
        dw    = z * gz;
        ndt   = numel(wave_exc);
        M_iso = zeros(135, N_iso, 4);

        % rotate pulse: +z -> +x
        wave_exc = wave_exc * exp(-1i * pi/2);
        
        % use mex based RK4 simulation
        parfor j = 1:135
            w1x  = real(wave_exc) / (abs(sum(wave_exc))*dt) * FAs(j);
            w1y  = imag(wave_exc) / (abs(sum(wave_exc))*dt) * FAs(j);
            for k=1:N_iso      
                M_iso(j,k,:) = mex_BLOCH_rk4([0; 0; 1; 1], w1x, w1y, dw(k) * ones(ndt, 1), 0, 0, dt);
            end
        end
        
        % split magnetization components
        Mxy = squeeze(M_iso(:,:,1)) + 1i * squeeze(M_iso(:,:,2));
        Mz  = squeeze(M_iso(:,:,3));
        
        % rephasing
        area_reph = -gz * ndt*dt/2;
        phi_reph  = repmat(z' * area_reph, [135, 1]);
        Mxy       = Mxy .* exp(-1i*phi_reph);
        
        % calculate slice profile
        [azimuthal, alpha, ~] = cart2sph(real(Mxy), imag(Mxy), Mz);
        alpha = pi/2 - alpha;
        SPROF_exc = alpha ./ repmat(FAs, [1, N_iso]) .* exp(1i*azimuthal);
        SPROF_exc = SPROF_exc.';
        
        % save slice profile simulation
        save(SProf_path, 'SPROF_exc', 'wave_exc', 'z', 'gz', 'dt');
    
    end

end

%% ------------------------------------------------------------------- TO DO !!!
function SPROF_rfc = presim_slice_profile_rfc(ID, GZ, wave_rfc, z, dt)

    % check equal gradient strength
    if numel(unique(GZ(ID==6)))>1
        error('different gradient strength detected!');
    else
        gz = unique(GZ(ID==6));
    end
    N_iso = numel(z);
    
    % search for slice profile simulation backup
    SProf_hash = ['slice_prof_rfc_' ...
                  num2str(numel(wave_rfc)) 'us_' ...
                  num2str(gz/1e3) 'kHzm_' ...
                  num2str(N_iso) 'Isos_' ...
                  num2str(min(z)*1e3) '_' num2str(max(z)*1e3) 'mm_' ...
                  pulseq_get_wave_hash(round(abs(wave_rfc),3))];
        
    SProf_hash  = [strrep(SProf_hash, '.', ',') '.mat'];
    SProf_path  = pulseq_get_path('MRF_find_slice_profile');
    SProf_path  = [SProf_path SProf_hash];
    SProf_exist = isfile(SProf_path);
    
    if SProf_exist==1
        load(SProf_path);
    else
        % calculate rf waveforms: 1° ... 135°
        % FAs = (136:1:270)' *pi/180;
        
        % bloch simulation for isochromats        
        % dw    = z * gz;
        % ndt   = numel(wave_rfc);
        % M_iso = zeros(135, N_iso, 4);      
        
        % use mex based RK4 simulation ... to do
        % parfor j = 1:135
        %     w1x  = real(wave_rfc) / (abs(sum(wave_rfc))*dt) * FAs(j);
        %     w1y  = imag(wave_rfc) / (abs(sum(wave_rfc))*dt) * FAs(j);
        %     for k=1:N_iso      
        %         M_iso(j,k,:) = mex_BLOCH_rk4([0; 0; 1; 1], w1x, w1y, dw(k) * ones(ndt, 1), 0, 0, dt);
        %     end
        % end           
    end

    SPROF_rfc = [];

    error('to do');

end
