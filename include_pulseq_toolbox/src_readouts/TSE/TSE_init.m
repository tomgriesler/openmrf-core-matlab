function [TSE, ktraj_adc, ktraj_full] = TSE_init(TSE, FOV, system, flag_plot)

% ---------------------------------------------------------
% ---------- init parameters and pulseq objects -----------
% ------------ readout: Turbo-Spin-Echo (TSE) -------------
% ---------------------------------------------------------

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% 0815 params
TSE.exc_flipangle = pi/2;
TSE.rfc_flipangle = pi;
TSE.exc_phase     = pi/2;
TSE.t_rise        = 250 *1e-6;
TSE.f_spR         = 1.0;
TSE.f_spS         = 0.5;
TSE.t_exwd        = TSE.exc_time + system.rfRingdownTime + system.rfDeadTime;
TSE.t_refwd       = TSE.rfc_time + system.rfRingdownTime + system.rfDeadTime;
TSE.t_Sp          = (TSE.TE - TSE.t_acq  - TSE.t_refwd) * 0.5;
TSE.t_Spex        = (TSE.TE - TSE.t_exwd - TSE.t_refwd) * 0.5;

%% slice excitation
if strcmp(TSE.mode_exc, 'sinc')
    [TSE.rf_exc, TSE.gz] = mr.makeSincPulse( TSE.exc_flipangle, ...
                                             system, ...
                                             'Duration', TSE.exc_time,...
                                             'SliceThickness', FOV.dz, ...
                                             'PhaseOffset', TSE.exc_phase, ...
                                             'timeBwProduct', TSE.exc_tbw, ...
                                             'apodization', 0.5, ...
                                             'use', 'excitation' );
elseif strcmp(TSE.mode_exc, 'sigpy_SLR')
    [TSE.rf_exc, TSE.gz] = SIGPY_SLR(TSE.exc_flipangle, TSE.exc_time, TSE.exc_phase, TSE.exc_tbw, 'ex', 'ls', 0.01, 0.01, 0 ,FOV.dz, system);
    TSE.rf_exc.use = 'excitation';
end

TSE.GSex = mr.makeTrapezoid( 'z', system, ...
                             'amplitude', TSE.gz.amplitude, ...
                             'FlatTime', TSE.t_exwd, ...
                             'riseTime', TSE.t_rise );
clear TSE.gz;

%% slice refocusing
if strcmp(TSE.mode_exc,'sinc')
    [TSE.rf_rfc, TSE.gz] = mr.makeSincPulse( TSE.rfc_flipangle, ...
                                             system, ...
                                             'Duration', TSE.rfc_time,...
                                             'SliceThickness', FOV.dz, ...
                                             'PhaseOffset', 0, ...
                                             'timeBwProduct', TSE.rfc_tbw, ...
                                             'apodization', 0.5, ...
                                             'use', 'refocusing' );
elseif strcmp(TSE.mode_rfc, 'sigpy_SLR')                            
    [TSE.rf_rfc, TSE.gz] = SIGPY_SLR(TSE.rfc_flipangle, TSE.rfc_time, 0, TSE.rfc_tbw, 'se', 'ls', 0.01, 0.01, 0 ,FOV.dz, system);
    TSE.rf_rfc.use = 'refocusing';
end                                

TSE.GSref = mr.makeTrapezoid( 'z', system, ...
                              'amplitude', TSE.gz.amplitude, ...
                              'FlatTime', TSE.t_refwd, ...
                              'riseTime', TSE.t_rise );
clear TSE.gz;

%%
TSE.AGSex  = TSE.GSex.area/2;
TSE.GSspr  = mr.makeTrapezoid( 'z', system, ...
                               'area', TSE.AGSex*(1+TSE.f_spS), ...
                               'duration', TSE.t_Sp, ...
                               'riseTime', TSE.t_rise );
                       
TSE.GSspex = mr.makeTrapezoid( 'z', system, ...
                               'area', TSE.AGSex*TSE.f_spS, ...
                               'duration', TSE.t_Spex, ...
                               'riseTime', TSE.t_rise );

%% frequency encoding
TSE.deltak_x = 1 / FOV.fov_x;
TSE.kWidth   = FOV.Nx * TSE.deltak_x;

TSE.GRacq = mr.makeTrapezoid( 'x', system, ...
                          'FlatArea', TSE.kWidth, ...
                          'FlatTime', TSE.t_acq, ...
                          'riseTime', TSE.t_rise );

% new git
if TSE.os_mode==0
    TSE.adcNSamples = FOV.Nx;
elseif TSE.os_mode==1
    TSE.adcNSamples = 2*FOV.Nx;
end
TSE.adcDwell    = round( (TSE.GRacq.flatTime-20e-6) / (TSE.adcNSamples) /100e-9 ) *100e-9;
TSE.adc         = mr.makeAdc( TSE.adcNSamples, 'Dwell', TSE.adcDwell);
TSE.adc.delay   = 10e-6;
                 
TSE.GRspr = mr.makeTrapezoid( 'x', system, ...
                          'area', TSE.GRacq.area * TSE.f_spR, ...
                          'duration', TSE.t_Sp, ...
                          'riseTime', TSE.t_rise );
                        
TSE.GRspex = mr.makeTrapezoid( 'x', system, ...
                           'area', TSE.GRacq.area * (1+TSE.f_spR), ...
                           'duration', TSE.t_Spex, ...
                           'riseTime', TSE.t_rise );

TSE.AGRspr   = TSE.GRspr.area;
TSE.AGRpreph = TSE.GRacq.area/2 + TSE.AGRspr;
TSE.GRpreph  = mr.makeTrapezoid( 'x', system, ...
                                 'Area', TSE.AGRpreph, ...
                                 'duration', TSE.t_Spex, ...
                                 'riseTime', TSE.t_rise );

%% phase encoding
TSE.deltak_y = 1/FOV.fov_y;

if strcmp(TSE.enc_mode, 'centric')
    [TSE.PEorder, TSE.nex, FOV.Ny] = TSE_CentricEncoding( FOV.Ny, TSE.n_echo );
end
if strcmp(TSE.enc_mode, 'linear')
    TSE.PEorder = reshape(-FOV.Ny/2:FOV.Ny/2-1, FOV.Ny/TSE.n_echo, TSE.n_echo)';
    TSE.nex     = size(TSE.PEorder,2);
end

TSE.phaseAreas = TSE.PEorder * TSE.deltak_y;

%% split gradients and recombine into blocks

% slice selection
TSE.GS1times = [0 TSE.GSex.riseTime];
TSE.GS1amp   = [0 TSE.GSex.amplitude];
TSE.GS1      = mr.makeExtendedTrapezoid( 'z', 'times', TSE.GS1times, 'amplitudes', TSE.GS1amp );

TSE.GS2times = [0 TSE.GSex.flatTime];
TSE.GS2amp   = [TSE.GSex.amplitude TSE.GSex.amplitude];
TSE.GS2      = mr.makeExtendedTrapezoid( 'z', 'times', TSE.GS2times, 'amplitudes', TSE.GS2amp );

TSE.GS3times = [0 TSE.GSspex.riseTime TSE.GSspex.riseTime + TSE.GSspex.flatTime TSE.GSspex.riseTime + TSE.GSspex.flatTime + TSE.GSspex.fallTime];
TSE.GS3amp   = [TSE.GSex.amplitude TSE.GSspex.amplitude TSE.GSspex.amplitude TSE.GSref.amplitude];
TSE.GS3      = mr.makeExtendedTrapezoid( 'z', 'times', TSE.GS3times, 'amplitudes', TSE.GS3amp );

TSE.GS4times = [0 TSE.GSref.flatTime];
TSE.GS4amp   = [TSE.GSref.amplitude TSE.GSref.amplitude];
TSE.GS4      = mr.makeExtendedTrapezoid( 'z', 'times', TSE.GS4times, 'amplitudes', TSE.GS4amp );

TSE.GS5times = [0 TSE.GSspr.riseTime TSE.GSspr.riseTime + TSE.GSspr.flatTime TSE.GSspr.riseTime + TSE.GSspr.flatTime + TSE.GSspr.fallTime];
TSE.GS5amp   = [TSE.GSref.amplitude TSE.GSspr.amplitude TSE.GSspr.amplitude 0];
TSE.GS5      = mr.makeExtendedTrapezoid( 'z', 'times', TSE.GS5times, 'amplitudes', TSE.GS5amp );

TSE.GS7times = [0 TSE.GSspr.riseTime TSE.GSspr.riseTime + TSE.GSspr.flatTime TSE.GSspr.riseTime + TSE.GSspr.flatTime + TSE.GSspr.fallTime];
TSE.GS7amp   = [0 TSE.GSspr.amplitude TSE.GSspr.amplitude TSE.GSref.amplitude];
TSE.GS7      = mr.makeExtendedTrapezoid( 'z', 'times', TSE.GS7times, 'amplitudes' , TSE.GS7amp );

% and now the readout gradient....

TSE.GR3 = TSE.GRpreph;%GRspex;

TSE.GR5times = [0 TSE.GRspr.riseTime TSE.GRspr.riseTime + TSE.GRspr.flatTime TSE.GRspr.riseTime + TSE.GRspr.flatTime + TSE.GRspr.fallTime];
TSE.GR5amp   = [0 TSE.GRspr.amplitude TSE.GRspr.amplitude TSE.GRacq.amplitude];
TSE.GR5      = mr.makeExtendedTrapezoid( 'x', 'times', TSE.GR5times, 'amplitudes', TSE.GR5amp );

TSE.GR6times = [0 TSE.t_acq];
TSE.GR6amp   = [TSE.GRacq.amplitude TSE.GRacq.amplitude];
TSE.GR6      = mr.makeExtendedTrapezoid( 'x', 'times', TSE.GR6times, 'amplitudes', TSE.GR6amp );

TSE.GR7times = [0 TSE.GRspr.riseTime TSE.GRspr.riseTime + TSE.GRspr.flatTime TSE.GRspr.riseTime + TSE.GRspr.flatTime + TSE.GRspr.fallTime];
TSE.GR7amp   = [TSE.GRacq.amplitude TSE.GRspr.amplitude TSE.GRspr.amplitude 0];
TSE.GR7      = mr.makeExtendedTrapezoid( 'x', 'times', TSE.GR7times, 'amplitudes', TSE.GR7amp );

% and filltimes
TSE.tex     = TSE.GS1.tt(end) + TSE.GS2.tt(end) + TSE.GS3.tt(end); % new git: .t->.tt
TSE.tref    = TSE.GS4.tt(end) + TSE.GS5.tt(end) + TSE.GS7.tt(end) + TSE.t_acq;  % new git: .t->.tt
TSE.tend    = TSE.GS4.tt(end) + TSE.GS5.tt(end); % new git: .t->.tt
TSE.tETrain = TSE.tex + TSE.n_echo * TSE.tref + TSE.tend;

%% rf: frequency and phase offsets
TSE.rfc_phase          = TSE.exc_phase - pi/2;
TSE.rf_exc.freqOffset  = TSE.GSex.amplitude   * FOV.z_offset;
TSE.rf_rfc.freqOffset  = TSE.GSref.amplitude  * FOV.z_offset;
TSE.rf_exc.phaseOffset = TSE.exc_phase - 2*pi * TSE.rf_exc.freqOffset * mr.calcRfCenter(TSE.rf_exc);
TSE.rfc_phase          = TSE.rfc_phase - 2*pi * TSE.rf_rfc.freqOffset * mr.calcRfCenter(TSE.rf_rfc);

%% prephaser
for loop_TR = 1 : TSE.nex
for loop_TE = 1 : TSE.n_echo
    TSE.GPpre(loop_TR,loop_TE).GPpre = mr.makeTrapezoid( 'y', system, 'Area', TSE.phaseAreas(loop_TE, loop_TR),  'Duration', TSE.t_Sp, 'riseTime', TSE.t_rise );
    TSE.GPrew(loop_TR,loop_TE).GPrew = mr.makeTrapezoid( 'y', system, 'Area', -TSE.phaseAreas(loop_TE, loop_TR), 'Duration', TSE.t_Sp, 'riseTime', TSE.t_rise );
end
end
clear loop TE loop_TR

%% TR filling delay
TSE.TR     = round(TSE.TR / system.gradRasterTime) * system.gradRasterTime;
TSE.t_objs = 0;
TSE.t_objs = TSE.t_objs + mr.calcDuration(TSE.GS1);
TSE.t_objs = TSE.t_objs + mr.calcDuration(TSE.GS2, TSE.rf_exc);
TSE.t_objs = TSE.t_objs + mr.calcDuration(TSE.GS3, TSE.GR3);
for loop_TE = 1 : TSE.n_echo
    TSE.t_objs = TSE.t_objs + mr.calcDuration(TSE.GS4, TSE.rf_rfc);
    TSE.t_objs = TSE.t_objs + mr.calcDuration(TSE.GS5, TSE.GR5, TSE.GPpre(1,loop_TE).GPpre);
    TSE.t_objs = TSE.t_objs + mr.calcDuration(TSE.GR6, TSE.adc);
    TSE.t_objs = TSE.t_objs + mr.calcDuration(TSE.GS7, TSE.GR7, TSE.GPrew(1,loop_TE).GPrew);    
end
TSE.t_objs = TSE.t_objs + mr.calcDuration(TSE.GS4);
TSE.t_objs = TSE.t_objs + mr.calcDuration(TSE.GS5);
clear loop TE

if TSE.t_objs >= TSE.TR-system.gradRasterTime
    TSE.tr_filling_delay = system.gradRasterTime;
    TSE.TR               = TSE.t_objs + system.gradRasterTime;
    disp(['TR increased to ' num2str(TSE.TR*1e3) 'ms !']);
else
    TSE.tr_filling_delay = TSE.TR - TSE.t_objs;
end

TSE.tr_filling_delay = round(TSE.tr_filling_delay / system.gradRasterTime) * system.gradRasterTime;
TSE.tr_filling_delay = mr.makeDelay(TSE.tr_filling_delay);

%% calculate kspace trajectory
seq = mr.Sequence(system);
for loop_TR = 1 : TSE.nex
    TSE_add();
end
[ktraj_adc, ktraj_full] = pulseq_get_ktraj(seq, flag_plot);

end

%% -----------------------------------------------------------------------------------
function [PEorder, nex, Ny] = TSE_CentricEncoding( Ny, necho )

    Nyold = Ny;
    
    nex = Ny / necho / 2;
    nex = ceil(nex) * 2;
    Ny  = nex * necho;
    dy  = Ny/2 / necho;
    dy  = ceil(dy);
    Ny  = dy * necho*2;
    
    
    pe_up   = (0:necho-1) * dy;
    pe_down = -pe_up - 1;
    
    PEorder = zeros(necho, nex);
    
    k = 1;
    for j=1:nex/2   
        PEorder(:,k)   = pe_up(:);
        PEorder(:,k+1) = pe_down(:);      
        pe_up          = pe_up+1;
        pe_down        = pe_down-1;
        k              = k+2;       
    end
    
    if Nyold ~= Ny
        display(['Ny changed!  ' num2str(Nyold) ' -> ' num2str(Ny)])
    end

end
