function [UTE, ktraj_adc, ktraj_full, ktraj_reco] = UTE_init(UTE, FOV, system)

% ---------------------------------------------------------
% ---------- init parameters and pulseq objects -----------
% ------------ readout: Ultra-short TE (UTE) --------------
% ---------------------------------------------------------

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
% source: https://pulseq.github.io/writeUTE_rs.html
% a basic UTE-like sequence
% achieves "TE" below 100 us

%% calculate projection angles
if strcmp(UTE.phi_mode, 'EqualPi')
    UTE.phi      = linspace(0, pi, UTE.NR+1);
    UTE.dphi     = UTE.phi(2) - UTE.phi(1);
    UTE.phi(end) = [];
end
if strcmp(UTE.phi_mode, 'Equal2Pi')
    UTE.phi      = linspace(0, 2*pi, UTE.NR+1);
    UTE.dphi     = UTE.phi(2) - UTE.phi(1);
    UTE.phi(end) = [];
end
if strcmp(UTE.phi_mode, 'RandomPi')
    UTE.phi       = linspace(0, pi, UTE.NR+1);
    UTE.dphi      = UTE.phi(2) - UTE.phi(1);    
    UTE.phi(end)  = [];
    [~, temp_ind] = sort(rand(1,UTE.NR));
    UTE.phi       = UTE.phi(temp_ind);
    clear temp_ind;
end
if strcmp(UTE.phi_mode, 'Random2Pi')
    UTE.phi       = linspace(0, 2*pi, UTE.NR+1);
    UTE.dphi      = UTE.phi(2) - UTE.phi(1);    
    UTE.phi(end)  = [];
    [~, temp_ind] = sort(rand(1,UTE.NR));
    UTE.phi       = UTE.phi(temp_ind);
    clear temp_ind;
end
if strcmp(UTE.phi_mode, 'GoldenAngle')
    if ~isfield(UTE, 'dphi')
        UTE.dphi = 2*pi / (1+sqrt(5));
    end    
    UTE.phi = (0:UTE.NR-1) * UTE.dphi;
end
if strcmp(UTE.phi_mode, 'RandomGoldenAngle')
    if ~isfield(UTE, 'dphi')
        UTE.dphi = 2*pi / (1+sqrt(5));
    end   
    UTE.phi       = (0:UTE.NR-1) * UTE.dphi;
    [~, temp_ind] = sort(rand(1,UTE.NR));
    UTE.phi       = UTE.phi(temp_ind);
    clear temp_ind;
end

if size(UTE.phi,1)<size(UTE.phi,2)
    UTE.phi = UTE.phi';
end       
if UTE.Ndummy>0
    UTE.phi = [zeros(UTE.Ndummy,1); UTE.phi];
end

%% Create alpha-degree slice selection pulse and gradient
if strcmp(UTE.exc_mode, 'sinc')
[rf, gz] = mr.makeSincPulse( UTE.exc_fa, 'system', system, ...
                             'Duration', UTE.exc_time, ...
                             'SliceThickness', FOV.dz, ...
                             'apodization', 0.5,...
                             'timeBwProduct', UTE.exc_tbw, ...
							 'use', 'excitation', ...
                             'centerpos', 1);
 elseif strcmp(UTE.exc_mode, 'sigpy_SLR')
        [rf, gz] = SIGPY_SLR( UTE.exc_fa, UTE.exc_time, 0, UTE.exc_tbw, 'st', 'ls', 0.01, 0.01, 0 , FOV.dz, system);
end

%% calculate gradients and adc

% resample the RF pulse to the ramp
gza         = [0 1 1 0];
gzt         = cumsum([0 gz.riseTime gz.flatTime gz.fallTime]);
gzas_0      = interp1(gzt + gz.delay, gza, rf.t + rf.delay);
rft_1       = [system.rfRasterTime:system.rfRasterTime:UTE.exc_time + 0.5 * gz.fallTime];
gzas_1      = interp1(gzt + gz.delay, gza, rft_1 + rf.delay + gz.fallTime * 0.5);
gzas_1(~isfinite(gzas_1)) = 0; % we are getting a NaN sometimes
kzs_0       = cumsum(gzas_0);
kzs_1       = cumsum(gzas_1);
kzs_0       = kzs_0 - max(kzs_0);
kzs_1       = kzs_1 - max(kzs_1);
rfs_1       = interp1(kzs_0, rf.signal, kzs_1);
rf.t        = rft_1;
rf.signal   = rfs_1.*gzas_1;
gz.flatTime = gz.flatTime - gz.fallTime * 0.5; % oops, we can get off gradient raster here, FIXME

% Define other gradients and ADC events
Nxo         = round(UTE.ro_os * FOV.Nxy);
deltak      = 1 / FOV.fov_xy / 2;
ro_area     = FOV.Nxy * deltak;
gx          = mr.makeTrapezoid('x', 'FlatArea', ro_area, 'FlatTime', UTE.adcTime, 'system', system);
adc_dur     = floor(gx.flatTime / Nxo * 1e7) * 1e-7 * Nxo; % round down dwell time to 100ns (Siemens ADC raster)
adc         = mr.makeAdc(Nxo, 'Duration', adc_dur, 'system', system);
gx.flatTime = gx.flatTime * UTE.spoil_ro;
UTE.TE      = ceil((UTE.TE + adc.dwell * UTE.ro_discard) / system.gradRasterTime) * system.gradRasterTime;
delayTR     = ceil((UTE.TR - mr.calcDuration(gz) - mr.calcDuration(gx) - UTE.TE) / system.gradRasterTime) * system.gradRasterTime;
gx.delay    = mr.calcDuration(gz) + UTE.TE;
adc.delay   = floor((gx.delay - adc.dwell * 0.5 - adc.dwell * UTE.ro_discard) / system.gradRasterTime) * system.gradRasterTime; % take into accout 0.5 samples ADC shift
fprintf('TE= %d us; delay in TR:= %d us\n', round(UTE.TE * 1e6), floor(delayTR * 1e6));

%% save objects
UTE.adc     = adc;
UTE.rf      = rf;
UTE.gx      = gx;
UTE.gz_pos  = gz;
UTE.gz_neg  = gz;
UTE.delayTR = mr.makeDelay(delayTR);
UTE.gz_neg.amplitude = -UTE.gz_neg.amplitude;
UTE.adcNSamples = Nxo;

%% save trajectory
seq = mr.Sequence(system);
for loop_NR = (1-UTE.Ndummy) : UTE.NR
    UTE_add();
end
[ktraj_adc, ktraj_full] = pulseq_get_ktraj(seq, 1);
 
ktraj_adc(3,:)    = [];
ktraj_full(3,:)   = [];
ktraj_x           = ktraj_adc(1,:);
ktraj_y           = ktraj_adc(2,:);
ktraj_x           = reshape(ktraj_x, UTE.adcNSamples, numel(ktraj_x)/UTE.adcNSamples)';
ktraj_y           = reshape(ktraj_y, UTE.adcNSamples, numel(ktraj_y)/UTE.adcNSamples)';
ktraj_reco(1,:,:) = ktraj_x(:,:);
ktraj_reco(2,:,:) = ktraj_y(:,:);

end
