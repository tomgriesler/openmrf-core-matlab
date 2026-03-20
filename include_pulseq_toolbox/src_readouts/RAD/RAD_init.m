function [RAD, ktraj_adc, ktraj_full, ktraj_reco] = RAD_init(RAD, FOV, system)

% ---------------------------------------------------------
% ---------- init parameters and pulseq objects -----------
% ----------------- readout: radial (RAD) -----------------
% ---------------------------------------------------------

% source: https://pulseq.github.io/writeFastRadialGradientEcho.html
% TE & TR are as short as possible
% version 20231112
% M. Gram: added rewinder gradients
%        corrected gxPrephaser -> check: Siemens sampling in the center of
%                                        the dwell periods!

%% calculate projection angles
if strcmp(RAD.phi_mode, 'EqualPi')
    RAD.phi      = linspace(0, pi, RAD.NR+1);
    RAD.dphi     = RAD.phi(2) - RAD.phi(1);
    RAD.phi(end) = [];
end
if strcmp(RAD.phi_mode, 'Equal2Pi')
    RAD.phi      = linspace(0, 2*pi, RAD.NR+1);
    RAD.dphi     = RAD.phi(2) - RAD.phi(1);
    RAD.phi(end) = [];
end
if strcmp(RAD.phi_mode, 'RandomPi')
    RAD.phi       = linspace(0, pi, RAD.NR+1);
    RAD.dphi      = RAD.phi(2) - RAD.phi(1);    
    RAD.phi(end)  = [];
    [~, temp_ind] = sort(rand(1,RAD.NR));
    RAD.phi       = RAD.phi(temp_ind);
    clear temp_ind;
end
if strcmp(RAD.phi_mode, 'Random2Pi')
    RAD.phi       = linspace(0, 2*pi, RAD.NR+1);
    RAD.dphi      = RAD.phi(2) - RAD.phi(1);    
    RAD.phi(end)  = [];
    [~, temp_ind] = sort(rand(1,RAD.NR));
    RAD.phi       = RAD.phi(temp_ind);
    clear temp_ind;
end
if strcmp(RAD.phi_mode, 'GoldenAngle')
    if ~isfield(RAD, 'dphi')
        RAD.dphi = 2*pi / (1+sqrt(5));
    end    
    RAD.phi = (0:RAD.NR-1) * RAD.dphi;
end
if strcmp(RAD.phi_mode, 'RandomGoldenAngle')
    if ~isfield(RAD, 'dphi')
        RAD.dphi = 2*pi / (1+sqrt(5));
    end   
    RAD.phi       = (0:RAD.NR-1) * RAD.dphi;
    [~, temp_ind] = sort(rand(1,RAD.NR));
    RAD.phi       = RAD.phi(temp_ind);
    clear temp_ind;
end

if size(RAD.phi,1)<size(RAD.phi,2)
    RAD.phi = RAD.phi';
end

%% calc adc
[RAD.adc, RAD.adcTime, RAD.adcNSamples, RAD.adcBW] = RAD_calc_adc(RAD.adcTime, FOV.Nxy*RAD.ro_os, system);

%% Create alpha-degree slice selection pulse and gradient
if strcmp(RAD.exc_mode, 'sinc')
[rf, gz, gzReph] = mr.makeSincPulse( RAD.exc_fa, 'system', system, ...
                                     'Duration', RAD.exc_time, ...
                                     'SliceThickness', FOV.dz, ...
                                     'apodization', 0.5,...
                                     'timeBwProduct', RAD.exc_tbw, ...
                                     'use', 'excitation' );
 elseif strcmp(RAD.exc_mode, 'sigpy_SLR')
        [rf, gz, gzReph] = SIGPY_SLR( RAD.exc_fa, RAD.exc_time, 0, RAD.exc_tbw, 'st', 'ls', 0.01, 0.01, 0 , FOV.dz, system);
end
gzReph.delay = mr.calcDuration(gz);
gzComb       = mr.addGradients({gz, gzReph}, 'system', system);

%% calculate readout gradient and read prephaser
deltak = 1/FOV.fov_xy;
gx     = mr.makeTrapezoid('x', 'FlatArea', FOV.Nxy*deltak, 'FlatTime', RAD.adcTime, 'system', system);

% calculate prephaser gradient (1st iteration)
gxPre = mr.makeTrapezoid('x', 'Area', -FOV.Nxy*deltak/2 - (gx.area-gx.flatArea)/2, 'system', system);
gxPre = mr.align('right', gxPre, 'right', gzComb);
gxPre = gxPre{1};
RAD.adc.delay = gx.riseTime;

% check kspace center position
seq = mr.Sequence(system);
seq.addBlock(rf, gzComb, gxPre);
seq.addBlock(gx, RAD.adc);
kx_corr = seq.calculateKspacePP();
kx_corr = kx_corr(1,RAD.adcNSamples/2+1);
%kx_corr = kx_corr + ?; Siemens sampling in the center of the dwell periods!
clear seq gxPre;

% calculate prephaser gradient (2nd iteration)
gxPre = mr.makeTrapezoid('x', 'Area', -FOV.Nxy*deltak/2 - (gx.area-gx.flatArea)/2 - kx_corr, 'system', system);
gxPre = mr.align('right', gxPre, 'right', gzComb);
gxPre = gxPre{1};
RAD.adc.delay = gx.riseTime;

% prevent gx dephasing during rf
addDelay = mr.calcDuration(rf) - gxPre.delay;
if addDelay > 0
    gxPre.delay = gxPre.delay+ceil(addDelay/system.gradRasterTime)*system.gradRasterTime;
end

%% gradient slice spoiling
if RAD.spoil_sl > 0
    sp_area_needed      = RAD.spoil_sl / FOV.dz - gz.area / 2;
    gzSpoil             = mr.makeTrapezoid('z', 'Area', sp_area_needed, 'system', system, 'Delay', gx.riseTime+gx.flatTime+system.gradRasterTime);
else
    gzSpoil = [];
end

%% rf spoiling
if strcmp(RAD.spoil_rf_mode, 'lin')
    RAD.spoil_rf_pow = 1;
elseif strcmp(RAD.spoil_rf_mode, 'quad')
    RAD.spoil_rf_pow = 2;
else
    error('spoiling increment mode unknown!');
end

%% gradient readout spoiling
if strcmp(RAD.mode_rewind, 'on')
    RAD.spoil_ro = 0;
end
if RAD.spoil_ro > 0
    ro_add_time = ceil(((gx.area/FOV.Nxy*(FOV.Nxy/2+1)*(1+RAD.spoil_ro))/gx.amplitude)/system.gradRasterTime)*system.gradRasterTime;
    gx.flatTime = gx.flatTime+ro_add_time; % careful, areas stored in the object are now wrong
end

%% calculate rewinder gradient
if strcmp(RAD.mode_rewind, 'on')
    seq = mr.Sequence(system);
    seq.addBlock(rf, gzComb, gxPre);
    seq.addBlock(gx, gzSpoil);
    [~, ~, temp_ktraj] = seq.calculateKspacePP();
    kx_final  = temp_ktraj(1,end);
    gx_rewind = mr.makeTrapezoid('x', 'Area', -kx_final, 'maxGrad', system.maxGrad*0.75, 'maxSlew', system.maxSlew*0.75, 'system', system);
    gx_rewind.delay = mr.calcDuration(gx);
    gxComb = mr.addGradients({gx, gx_rewind}, 'system', system);
    clear seq;
else
    gxComb = gx;
end

%% save RAD objects
RAD.rf      = rf;
RAD.gzComb  = gzComb;
RAD.gxPre   = gxPre;
RAD.gxComb  = gxComb;
RAD.gzSpoil = gzSpoil;

%% calc TR filling delay
temp_tr = mr.calcDuration(RAD.rf, RAD.gzComb, RAD.gxPre) + ...
          mr.calcDuration(RAD.gxComb, RAD.gzSpoil, RAD.adc);
if temp_tr > RAD.TR
    RAD.TR = temp_tr + system.gradRasterTime;
    warning(['TR corrected to ' num2str(RAD.TR*1e3) 'ms']);
    RAD.tr_filling_delay = system.gradRasterTime;
else
    RAD.tr_filling_delay = RAD.TR - temp_tr;
end
RAD.tr_filling_delay = round(RAD.tr_filling_delay/system.gradRasterTime) * system.gradRasterTime;
RAD.tr_filling_delay = mr.makeDelay(RAD.tr_filling_delay);

%% save trajectory
seq = mr.Sequence(system);
for loop_NR = (1-RAD.Ndummy) : RAD.NR
    RAD_add();
end
[ktraj_adc, ktraj_full] = pulseq_get_ktraj(seq, 1);
 
ktraj_adc(3,:)    = [];
ktraj_full(3,:)   = [];
ktraj_x           = ktraj_adc(1,:);
ktraj_y           = ktraj_adc(2,:);
ktraj_x           = reshape(ktraj_x, RAD.adcNSamples, numel(ktraj_x)/RAD.adcNSamples)';
ktraj_y           = reshape(ktraj_y, RAD.adcNSamples, numel(ktraj_y)/RAD.adcNSamples)';
ktraj_reco(1,:,:) = ktraj_x(:,:);
ktraj_reco(2,:,:) = ktraj_y(:,:);

end
