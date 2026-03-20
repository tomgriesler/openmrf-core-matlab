function [GRE, ktraj_adc, ktraj_full] = GRE_init(GRE, FOV, system)

% ---------------------------------------------------------
% ---------- init parameters and pulseq objects -----------
% ------------- readout: GRadient-Echo (GRE) --------------
% ---------------------------------------------------------

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% limit max grad strength and slew rate
if ~isfield(GRE, 'max_grad')
    GRE.max_grad   = 0.95;
end
if ~isfield(GRE, 'max_slew')
    GRE.max_slew   = 0.75;
end
system.maxGrad = system.maxGrad * GRE.max_grad;
system.maxSlew = system.maxSlew * GRE.max_slew;

%% check number of different FAs, TEs and TRs
GRE.n_FAs = numel(GRE.exc_flipangle);  % different FAs for Fingerprinting
GRE.n_TRs = numel(GRE.TRs);            % different TRs for Fingerprinting
GRE.n_TEs = numel(GRE.TEs);            % different TEs for B0 mapping

% use colum vectors
if size(GRE.exc_flipangle,2) > size(GRE.exc_flipangle,1)
    GRE.exc_flipangle = GRE.exc_flipangle';
end
if size(GRE.TEs,2) > size(GRE.TEs,1)
    GRE.TEs = GRE.TEs';
end
if size(GRE.TRs,2) > size(GRE.TRs,1)
    GRE.TRs = GRE.TRs';
end

% check list size
if GRE.n_TRs > 1  &&  GRE.n_TEs == 1
    GRE.TEs = GRE.TEs * ones(GRE.n_TRs, 1);
end
if GRE.n_TEs > 1  &&  GRE.n_TRs == 1
    GRE.TRs = GRE.TRs * ones(GRE.n_TEs, 1);
end
if GRE.n_TRs > 1  &&  GRE.n_TEs > 1  &&  GRE.n_TRs ~= GRE.n_TEs
    error('use identical number of TEs and TRs !');
end

GRE.n_TRs = numel(GRE.TRs);
GRE.n_TEs = numel(GRE.TEs);

%% calculate rf objects
for j = 1 : GRE.n_FAs
    if strcmp(GRE.exc_mode, 'sinc')
        [GRE.rf(j), GRE.gz] = mr.makeSincPulse( GRE.exc_flipangle(j), ...
                                                system, ...
                                                'Duration', GRE.exc_time, ...
                                                'timeBwProduct', GRE.exc_tbw, ...
                                                'apodization', 0.5, ...
                                                'SliceThickness', FOV.dz, ...
                                                'use', 'excitation'); 
        GRE.exc_shape = [];
    elseif strcmp(GRE.exc_mode, 'sigpy_SLR')
        [GRE.rf(j), GRE.gz] = SIGPY_SLR( GRE.exc_flipangle(j), GRE.exc_time, 0, GRE.exc_tbw, GRE.exc_shape, 'ls', 0.01, 0.01, 0 , FOV.dz, system);
    end
    GRE.rf(j).freqOffset  = GRE.gz.amplitude * FOV.z_offset;
    GRE.rf(j).phaseOffset = 0; % changed by phase cycling
end
clear j;

%% k-space parameters
GRE.deltak_x   = 1/FOV.fov_x;
GRE.deltak_y   = 1/FOV.fov_y;
GRE.gx         = mr.makeTrapezoid('x', 'FlatArea', FOV.Nx*GRE.deltak_x, 'FlatTime', GRE.t_acq, 'system', system);
GRE.phaseAreas = ( -FOV.Ny/2 : FOV.Ny/2-1 ) * GRE.deltak_y;

%% read and phase prephaser, slice rephaser
if isempty(GRE.t_pre) % auto-minimization of t_pre
    gdx = mr.makeTrapezoid('x', 'Area', -GRE.gx.area/2,           'maxGrad', system.maxGrad, 'maxSlew', system.maxSlew, 'system', system);
    gdy = mr.makeTrapezoid('y', 'Area', max(abs(GRE.phaseAreas)), 'maxGrad', system.maxGrad, 'maxSlew', system.maxSlew, 'system', system);
    gdz = mr.makeTrapezoid('z', 'Area', -GRE.gz.area/2,           'maxGrad', system.maxGrad, 'maxSlew', system.maxSlew, 'system', system);
    GRE.t_pre = max([ mr.calcDuration(gdx), mr.calcDuration(gdy), mr.calcDuration(gdz) ]);
    clear gdx gdy gdz;
end
GRE.gx_pre = mr.makeTrapezoid('x', 'Area', -GRE.gx.area/2, 'Duration', GRE.t_pre, 'system', system);
for j = 1 : FOV.Ny
    GRE.gy_pre(j) = mr.makeTrapezoid('y', 'Area', GRE.phaseAreas(j), 'Duration', GRE.t_pre, 'system', system);
end
GRE.gz_reph = mr.makeTrapezoid('z', 'Area', -GRE.gz.area/2, 'Duration', GRE.t_pre, 'system', system);

%% read and phase rewinder, slice spoiler or rewinder

% check slice spoiler mode
if abs(GRE.spoil_nTwist) > 0
    GRE.spoil_gz_mode = 'unbalanced';
else
    GRE.spoil_gz_mode = 'balanced';
end

 % auto-minimization of t_spoil
if isempty(GRE.t_spoil)
    gdx = mr.makeTrapezoid('x', 'Area', -GRE.gx.area/2,           'maxGrad', system.maxGrad, 'maxSlew', system.maxSlew, 'system', system);
    gdy = mr.makeTrapezoid('y', 'Area', max(abs(GRE.phaseAreas)), 'maxGrad', system.maxGrad, 'maxSlew', system.maxSlew, 'system', system);
    if strcmp(GRE.spoil_gz_mode, 'balanced')
        gdz = mr.makeTrapezoid('z', 'Area', -GRE.gz.area/2, 'maxGrad', system.maxGrad, 'maxSlew', system.maxSlew, 'system', system);
    elseif strcmp(GRE.spoil_gz_mode, 'unbalanced')
        gdz = mr.makeTrapezoid('z', 'Area', GRE.spoil_nTwist/FOV.dz - GRE.gz.area/2, 'maxGrad', system.maxGrad, 'maxSlew', system.maxSlew, 'system', system);
    end
    GRE.t_spoil = max([ mr.calcDuration(gdx), mr.calcDuration(gdy), mr.calcDuration(gdz) ]);
    clear gdx gdy gdz;
end

% calc x y rewinder
GRE.gx_rew = mr.makeTrapezoid('x', 'Area', -GRE.gx.area/2, 'Duration', GRE.t_spoil, 'system', system);
for j = 1 : FOV.Ny
    GRE.gy_rew(j) = mr.makeTrapezoid('y', 'Area', -GRE.phaseAreas(j), 'Duration', GRE.t_spoil, 'system', system);
end

%% slice spoiler
if abs(GRE.spoil_nTwist) > 0
    GRE.gz_spoil = mr.makeTrapezoid('z', 'Area', GRE.spoil_nTwist/FOV.dz - GRE.gz.area/2, 'Duration', GRE.t_spoil, 'system', system);
else
    GRE.gz_spoil = mr.makeTrapezoid('z', 'Area', -GRE.gz.area/2,                          'Duration', GRE.t_spoil, 'system', system);    
end

%% rf spoiling
if strcmp(GRE.spoil_rf_mode, 'quad')
    GRE.spoil_rf_pow = 2;
elseif strcmp(GRE.spoil_rf_mode, 'lin')
    GRE.spoil_rf_pow = 1;
else
    error('unknown rf spoiling mode!')
end

%% calculate adc
if GRE.os_mode==0
    GRE.adcNSamples = FOV.Nx;
elseif GRE.os_mode==1
    GRE.adcNSamples = 2*FOV.Nx;
end
GRE.adcDwell  = round( GRE.gx.flatTime / GRE.adcNSamples /100e-9 ) *100e-9;
GRE.adcBW     = 1 / GRE.adcDwell;
GRE.adc       = mr.makeAdc(GRE.adcNSamples, 'Dwell', GRE.adcDwell); 
GRE.adc.delay = round(GRE.gx.riseTime/system.gradRasterTime) * system.gradRasterTime;

%% calculate timings: filling delays
changed_TE_TR = 0;

for j=1:GRE.n_TEs

% try TE filling delay
GRE.delayTE(j) = GRE.TEs(j) ...
               - GRE.gz.flatTime/2 ...
               - GRE.gz.fallTime ...
               - mr.calcDuration(GRE.gx_pre, GRE.gy_pre(1), GRE.gz_reph) ...
               - mr.calcDuration(GRE.gx)/2;

% check TE filling delay
if GRE.delayTE(j) <= 0
    GRE.TEs(j) = 0;
end

% minimize TE filling delay
if GRE.TEs(j) == 0
    GRE.delayTE(j) = system.gradRasterTime;
    GRE.TEs(j)     = GRE.gz.flatTime/2 ...
                   + GRE.gz.fallTime ...
                   + mr.calcDuration(GRE.gx_pre, GRE.gy_pre(1), GRE.gz_reph) ...
                   + mr.calcDuration(GRE.gx)/2 ...
                   + system.gradRasterTime;
    changed_TE_TR  = 1;
end
GRE.delayTE(j) = ceil(GRE.delayTE(j) / system.gradRasterTime) * system.gradRasterTime;

% try TR filling delay
GRE.delayTR(j) = GRE.TRs(j) ...
               - mr.calcDuration(GRE.gz) ...
               - mr.calcDuration(GRE.gx_pre, GRE.gy_pre(1), GRE.gz_reph) ...
               - mr.calcDuration(GRE.gx) ...
               - mr.calcDuration(GRE.gx_rew, GRE.gy_rew(1), GRE.gz_spoil) ...
               - GRE.delayTE(j);

% check TR filling delay
if GRE.delayTR(j) <= 0
    GRE.TRs(j) = 0;
end

% minimize TR filling delay
if GRE.TRs(j) == 0
    GRE.delayTR(j) = system.gradRasterTime;
    GRE.TRs(j)     = mr.calcDuration(GRE.gz) ... 
                   + mr.calcDuration(GRE.gx_pre, GRE.gy_pre(1), GRE.gz_reph) ...
                   + mr.calcDuration(GRE.gx) ...
                   + mr.calcDuration(GRE.gx_rew, GRE.gy_rew(1), GRE.gz_spoil) ...
                   + max(GRE.delayTE(j)) ...
                   + system.gradRasterTime;
    changed_TE_TR  = 1;
end
GRE.delayTR(j) = ceil(GRE.delayTR(j) / system.gradRasterTime) * system.gradRasterTime;

end

if changed_TE_TR==1
    disp(['changed TEs to: [ms]']);
    disp(round(GRE.TEs*1e3,2));
    disp(['changed TRs to: [ms]']);
    disp(round(GRE.TRs*1e3,2));
end

%% calculate kspace trajectory
seq = mr.Sequence(system);
for loop_Ny = 1:FOV.Ny
    GRE_add();
end
[ktraj_adc, ktraj_full] = pulseq_get_ktraj(seq, 1);

end
