function SL = SL_init(SL, FOV, system)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% function for spin-lock initialization
% - create descriptor for spin-lock preparation
% - set default values
% - check tSL and fSL values
% - check remaining parameters
% - prevent timing errors
% - init SL objects
% - create main SL objects

    %% create descriptor for spin-lock preparation
    SL.info = {
        'seq_type:',       '[string] spin lock module type, >SSL<, >hSSL<, >qSSL<, >RESL<, >CSL<, >PSCSL<, >TBSL< or >BSL<';
        'relax_type',      '[string] relaxation time mode: >T1p<, >T2p<, >T2<';
        'tSL:',            '[s]      spin lock time array';
        'fSL:',            '[Hz]     spin lock amplitude array';
        'inv_mode:',       '[string] Mz or -Mz preparation >on< or >off<';
        'rho_axis:',       '[rad]    direction of locking field';
        'exc_mode:',       '[string] excitation mode: >sinc<, >adiabatic_AHP<, >adiabatic_BIR4<, >bp>, or >sigpy_SLR<';
        'exc_time:',       '[s]      excitation time';
        'exc_flip_angle:', '[rad]    should be pi/2';
        'rfc_mode:',       '[string] refocusing mode: >sinc<, >bp<, >comp< or >sigpy_SLR<';
        'rfc_time:',       '[s]      refocusing time';
        'rfc_flip_angle:', '[rad]    should be pi';
        'adia_wmax:',      '[rad/s]  adiabatic maximum w1';
    };
    
    %% set default values
    if ~isfield(SL, 'relax_type')
        SL.relax_type = 'T1p';
    end
    if ~isfield(SL, 'seq_type')
        SL.seq_type = 'BSL';
    end
    if ~isfield(SL, 'rho_axis')
        SL.rho_axis = 0;
    end
    if ~isfield(SL, 'freqOffset')
        SL.freqOffset = 0;
    end
    if ~isfield(SL, 'inv_mode')
        SL.inv_mode = {'off'};
    end
    if ~isfield(SL, 'exc_mode')
        SL.exc_mode = 'adiabatic_AHP';
    end
    if strcmp(SL.exc_mode, 'adiabatic_AHP')
        if ~isfield(SL, 'exc_time')
            SL.exc_time = 3.0 * 1e-3;
        end
        if ~isfield(SL, 'adia_wmax')
            SL.adia_wmax = 600 * 2*pi;    
        end
    end
    if strcmp(SL.exc_mode, 'sinc') || strcmp(SL.exc_mode, 'sigpy_SLR')
        if ~isfield(SL, 'exc_flip_angle')
            SL.exc_flip_angle = pi/2;
        end
        if ~isfield(SL, 'exc_time')
            SL.exc_time = 2.0 * 1e-3;
        end
        if ~isfield(SL, 'exc_tbw')
            SL.exc_tbw = 4;    
        end
    end
    if strcmp(SL.exc_mode, 'bp')
        if ~isfield(SL, 'exc_flip_angle')
            SL.exc_flip_angle = pi/2;
        end
        if ~isfield(SL, 'exc_time')
            SL.exc_time = 1.0 * 1e-3;
        end
    end
    if ~isfield(SL, 'rfc_mode')
        SL.rfc_mode = 'bp';
    end
    if strcmp(SL.rfc_mode, 'bp')
        if ~isfield(SL, 'rfc_time')
            SL.rfc_time = 1.0 * 1e-3;
        end
        if ~isfield(SL, 'rfc_flip_angle')
            SL.rfc_flip_angle = pi;
        end
    end
    if strcmp(SL.rfc_mode, 'sinc') || strcmp(SL.rfc_mode, 'sigpy_SLR')
        if ~isfield(SL, 'rfc_time')
            SL.rfc_time = 4.0 * 1e-3;
        end
        if ~isfield(SL, 'rfc_flip_angle')
            SL.rfc_flip_angle = pi;
        end
        if ~isfield(SL, 'rfc_tbw')
            SL.rfc_tbw = 4;
        end
    end
    if strcmp(SL.rfc_mode, 'comp')
        if ~isfield(SL, 'rfc_time')
            SL.rfc_time = 2.0 * 1e-3;
        end
    end
    if ~isfield(SL, 'crush_nTwists_x')
        SL.crush_nTwists_x = 4;   % [] number of 2pi twists in x direction
    end
    if ~isfield(SL, 'crush_nTwists_y')
        SL.crush_nTwists_y = 4;   % [] number of 2pi twists in y direction
    end
    if ~isfield(SL, 'crush_nTwists_z')
        SL.crush_nTwists_z = 11.3;   % [] number of 2pi twists in z direction
    end
    if ~isfield(SL, 'crush_lim_grad')
        SL.crush_lim_grad = 1/sqrt(3); % reduce crusher gradient amplitude from nominal limit
    end    
    if ~isfield(SL, 'crush_lim_slew')
        SL.crush_lim_slew = 1/sqrt(3); % reduce crusher gradient slew rate from nominal limit to avoid stimulation
    end

    %% clear unnecessary parameters
    if strcmp(SL.exc_mode, 'adiabatic_AHP')
        SL.exc_flip_angle = [];
        SL.exc_tbw = [];
    end
    if strcmp(SL.exc_mode, 'sinc') || strcmp(SL.exc_mode, 'sigpy_SLR') || strcmp(SL.exc_mode, 'bp')
        SL.adia_wmax = [];
    end
    if strcmp(SL.exc_mode, 'bp')
        SL.exc_tbw = [];
    end
    if strcmp(SL.rfc_mode, 'comp')
        SL.rfc_flip_angle = [];
        SL.rfc_tbw = [];
    end
    if strcmp(SL.rfc_mode, 'bp')
        SL.rfc_tbw = [];
    end

    %% check tSL and fSL values
    if size(SL.tSL,1) < size(SL.tSL,2)
        SL.tSL = SL.tSL';  % use column vectors
    end
    if size(SL.fSL,1) < size(SL.fSL,2)
        SL.fSL = SL.fSL';  % use column vectors
    end
    if numel(SL.tSL) ~= numel(SL.fSL)
        error('use same number of tSL and fSL values!')
    end
    if min(SL.tSL) < system.gradRasterTime * 4
        error('increase tSL!')
    end
    if max(SL.tSL) > 400 *1e-3
        warning('tSL too large!')
    end
    if max(SL.fSL) > 2000
        warning('fSL too large!')
    end
       
    %% prevent timing errors
    SL.tSL      = round(SL.tSL/system.rfRasterTime/4)    * system.rfRasterTime * 4;
    SL.exc_time = round(SL.exc_time/system.rfRasterTime) * system.rfRasterTime;
    SL.rfc_time = round(SL.rfc_time/system.rfRasterTime) * system.rfRasterTime;
    
    %% init SL objects
    
    SL.nSL = numel(SL.tSL); % number of different SL preparations
    
    % check list of relax types
    if numel(SL.relax_type)==1
        for j=1:SL.nSL
            SL.relax_type{j,1} = SL.relax_type{1};
        end
    end
    if size(SL.relax_type,1) < size(SL.relax_type,2)
        SL.relax_type = SL.relax_type';
    end
    if numel(SL.relax_type) ~= SL.nSL
        error('wrong number of relax types!');
    end
    
    % check list of seq types
    if numel(SL.seq_type)==1
        for j=1:SL.nSL
            SL.seq_type{j,1} = SL.seq_type{1};
        end
    end
    if size(SL.seq_type,1) < size(SL.seq_type,2)
        SL.seq_type = SL.seq_type';
    end
    if numel(SL.seq_type) ~= SL.nSL
        error('wrong number of seq types!');
    end
    
    % check list of inversion modes
    if numel(SL.inv_mode)==1
        for j=1:SL.nSL
            SL.inv_mode{j,1} = SL.inv_mode{1};
        end
    end
    if size(SL.inv_mode,1) < size(SL.inv_mode,2)
        SL.inv_mode = SL.inv_mode';
    end
    if numel(SL.inv_mode) ~= SL.nSL
        error('wrong number of inversion modes!');
    end
    
    % init SL main objects
    for j = 1:SL.nSL
        SL.SL_objs(j).relax_type = SL.relax_type{j};
        SL.SL_objs(j).seq_type   = SL.seq_type{j};
        SL.SL_objs(j).tSL        = SL.tSL(j);
        SL.SL_objs(j).fSL        = SL.fSL(j);
        SL.SL_objs(j).inv_mode   = SL.inv_mode{j};
        SL.SL_objs(j).exc_mode   = SL.exc_mode;
    end
    
    %% create main SL objects
    for j = 1:SL.nSL
        [SL.SL_objs(j).EXC1, SL.SL_objs(j).EXC2]                                 = SL_get_objs_EXC(SL, system, j);
        [SL.SL_objs(j).RFC1, SL.SL_objs(j).RFC2]                                 = SL_get_objs_RFC(SL, system, j);
        [SL.SL_objs(j).SL1, SL.SL_objs(j).SL2, SL.SL_objs(j).SL3]                = SL_get_objs_SL(SL, system, j);
        [SL.SL_objs(j).gx_crush, SL.SL_objs(j).gy_crush, SL.SL_objs(j).gz_crush] = CRUSH_x_y_z(SL.crush_nTwists_x, SL.crush_nTwists_y, SL.crush_nTwists_z, FOV.dx, FOV.dy, FOV.dz, SL.crush_lim_grad, SL.crush_lim_slew, system);
        SL.SL_objs(j).d1                                                         = mr.makeDelay(system.gradRasterTime);
        SL.SL_objs(j).d2                                                         = mr.makeDelay(system.gradRasterTime);
        SL.SL_objs(j).duration                                                   = SL_get_duration(SL.SL_objs(j));
    end

end

%% ----------------------------------------------------------------------------------------
function [EXC1, EXC2] = SL_get_objs_EXC(SL, system, loop_SL)

    % Global Excitation Pulses for SL Preparation
    
    % switch cases for pulse phases
    rho_axis = SL.rho_axis;
    
    if strcmp(SL.exc_mode, 'adiabatic_BIR4')
        error('to do!')
    end
    
    % T1p, T2, MZREX, OREX, SIRS case
    if ~strcmp(SL.exc_mode, 'adiabatic_AHP') && strcmp(SL.inv_mode(loop_SL), 'off')
        exc_phase1 = rho_axis - pi/2;
        exc_phase2 = rho_axis + pi/2;
    end
    if ~strcmp(SL.exc_mode, 'adiabatic_AHP') && strcmp(SL.inv_mode(loop_SL), 'on')
        exc_phase1 = rho_axis - pi/2;
        exc_phase2 = rho_axis - pi/2;
    end
    if strcmp(SL.exc_mode, 'adiabatic_AHP') && strcmp(SL.inv_mode(loop_SL), 'off')
        exc_phase1 = rho_axis;
        exc_phase2 = rho_axis;
    end
    if strcmp(SL.exc_mode, 'adiabatic_AHP') && strcmp(SL.inv_mode(loop_SL), 'on')
        exc_phase1 = rho_axis;
        exc_phase2 = rho_axis + pi;
    end
    
    % T2p case
    if strcmp(SL.relax_type(loop_SL), 'T2p')
    if ~strcmp(SL.exc_mode, 'adiabatic_AHP') && strcmp(SL.inv_mode(loop_SL), 'off')
        exc_phase1 = rho_axis;
        exc_phase2 = rho_axis + pi;
    end
    if ~strcmp(SL.exc_mode, 'adiabatic_AHP') && strcmp(SL.inv_mode(loop_SL), 'on')
        exc_phase1 = rho_axis;
        exc_phase2 = rho_axis;
    end
    if strcmp(SL.exc_mode, 'adiabatic_AHP') && strcmp(SL.inv_mode(loop_SL), 'off')
        exc_phase1 = rho_axis + pi/2;
        exc_phase2 = rho_axis + pi/2;
    end
    if strcmp(SL.exc_mode, 'adiabatic_AHP') && strcmp(SL.inv_mode(loop_SL), 'on')
        exc_phase1 = rho_axis + pi/2;
        exc_phase2 = rho_axis - pi/2;
    end
    end
    
    % calc pulses
    if strcmp(SL.exc_mode, 'sinc')
        EXC1 = mr.makeSincPulse( SL.exc_flip_angle, ...
                                 system, ...
                                 'Duration', SL.exc_time,...
                                 'timeBwProduct' , SL.exc_tbw, ...
                                 'SliceThickness', 1e6, ...
                                 'PhaseOffset', exc_phase1, ...
							     'use', 'excitation' ); 
        EXC2 = mr.makeSincPulse( SL.exc_flip_angle, ...
                                 system, ...
                                 'Duration', SL.exc_time,...
                                 'timeBwProduct' , SL.exc_tbw, ...
                                 'SliceThickness', 1e6, ...
                                 'PhaseOffset', exc_phase2, ...
							     'use', 'excitation' );
    elseif strcmp(SL.exc_mode, 'bp')
        EXC1 = SL_get_block_pulse(system, SL.exc_time, SL.exc_flip_angle, exc_phase1, 0, 'excitation');    
        EXC2 = SL_get_block_pulse(system, SL.exc_time, SL.exc_flip_angle, exc_phase2, 0, 'excitation');
    
    elseif strcmp(SL.exc_mode, 'adiabatic_AHP')
        AHP  = AHP_get_params(SL.exc_time, SL.adia_wmax);
        EXC1 = AHP_pulse( system, 'tipdown', AHP.tau, AHP.wmax, AHP.ratio, AHP.beta1, AHP.beta2,  exc_phase1 );
        EXC2 = AHP_pulse( system, 'tipup',   AHP.tau, AHP.wmax, AHP.ratio, AHP.beta1, AHP.beta2,  exc_phase2 );
    
    elseif strcmp(SL.exc_mode, 'sigpy_SLR')
        EXC1 = SIGPY_SLR( SL.exc_flip_angle, SL.exc_time, exc_phase1, SL.exc_tbw, 'ex', 'ls', 0.001, 0.001, 0 ,1e6, system);
        EXC2 = SIGPY_SLR( SL.exc_flip_angle, SL.exc_time, exc_phase2, SL.exc_tbw, 'ex', 'ls', 0.001, 0.001, 0 ,1e6, system);
    end

end

%% ----------------------------------------------------------------------------------------
function [RFC1, RFC2] = SL_get_objs_RFC(SL, system, loop_SL)

    % Global Refocusing Pulses for SL Preparation
    
    % calculate refocusing pulses
    rfc_phase1 = SL.rho_axis;
    rfc_phase2 = SL.rho_axis + pi;
    
    if strcmp(SL.rfc_mode, 'sinc')
        RFC1 = mr.makeSincPulse( SL.rfc_flip_angle, ...
                                 system, ...
                                 'Duration', SL.rfc_time,...
                                 'timeBwProduct', SL.rfc_tbw, ...
                                 'SliceThickness', 1e6, ...
                                 'PhaseOffset', rfc_phase1, ...
							     'use', 'refocusing' );    
        
        RFC2 = mr.makeSincPulse( SL.rfc_flip_angle, ...
                                 system, ...
                                 'Duration', SL.rfc_time,...
                                 'timeBwProduct', SL.rfc_tbw, ...
                                 'SliceThickness', 1e6, ...
                                 'PhaseOffset', rfc_phase2, ...
							     'use', 'refocusing' );
    
    elseif strcmp(SL.rfc_mode, 'bp')
        RFC1 = SL_get_block_pulse(system, SL.rfc_time, SL.rfc_flip_angle, rfc_phase1, 0, 'refocusing');    
        RFC2 = SL_get_block_pulse(system, SL.rfc_time, SL.rfc_flip_angle, rfc_phase2, 0, 'refocusing');
    
    elseif strcmp(SL.rfc_mode, 'sigpy_SLR')
        RFC1 = SIGPY_SLR( SL.rfc_flip_angle, SL.rfc_time, rfc_phase1, SL.rfc_tbw, 'se', 'ls', 0.001, 0.001, 0 ,1e6, system);
        RFC2 = SIGPY_SLR( SL.rfc_flip_angle, SL.rfc_time, rfc_phase2, SL.rfc_tbw, 'se', 'ls', 0.001, 0.001, 0 ,1e6, system);
    
    elseif strcmp(SL.rfc_mode, 'comp')
        RFC1 = SL_get_comp_pulse(system, SL.rfc_time, rfc_phase1);
        RFC2 = SL_get_comp_pulse(system, SL.rfc_time, rfc_phase2);
    
    end
    
    % switch cases for seq types
    if strcmp(SL.seq_type(loop_SL), 'SSL') || strcmp(SL.seq_type(loop_SL), 'hSSL') || strcmp(SL.seq_type(loop_SL), 'qSSL') || strcmp(SL.seq_type(loop_SL), 'RESL')
        RFC1 = [];
        RFC2 = [];
    end
    if strcmp(SL.seq_type(loop_SL), 'CSL') || strcmp(SL.seq_type(loop_SL), 'PSCSL')
        RFC2 = [];
    end

end

%% ----------------------------------------------------------------------------------------
function [SL1, SL2, SL3] = SL_get_objs_SL(SL, system, loop_SL)

    % ----- function for calculating different SL prep modules: -----
    % SSL:   Sepponen et al. 1985. doi.org/10.1097/00004728-198511000-00002
    % RESL:  Charagundla et al. 2003. doi.org/10.1016/S1090-7807(02)00197-0
    % CSL:   Witschey et al. 2007 doi.org/10.1016/j.jmr.2007.01.015
    % PSCSL: Bogdan et al. 2016. doi.org/10.1016/j.jmr.2016.04.017
    % TBSL:  Gram et al. 2019. https://archive.ismrm.org/2019/1215.html
    % BSL:   Gram et al. 2020. doi.org/10.1002/mrm.28585
    
    % note: hSSL (half-SSL) works like a SSL module, but the SL time
    % is subdivided in two parts. This helps for long SL times, 
    % qSSL works like hSSL with a subdivision in 4 parts
    
    % note: Witschey's CSL module needs the inversion option for the tip up pulse
    
    % calculate spin lock pulses
    SL1 = [];
    SL2 = [];
    SL3 = [];
    
    seq_type = SL.seq_type(loop_SL);
    tSL      = SL.tSL(loop_SL);
    fSL      = SL.fSL(loop_SL);
    phase1   = SL.rho_axis;
    phase2   = SL.rho_axis + pi;
    
    % use block pulses
    if tSL>50e-6
    if strcmp(seq_type, 'SSL')
        SL1 = SL_get_SL_pulse(system, tSL, fSL, phase1, SL.freqOffset);
    
    elseif strcmp(seq_type, 'hSSL')
        SL1 = SL_get_SL_pulse(system, tSL/2, fSL, phase1, SL.freqOffset);
        SL2 = SL1;
    
    elseif strcmp(seq_type, 'qSSL')
        SL1 = SL_get_SL_pulse(system, tSL/4, fSL, phase1, SL.freqOffset);
        SL2 = SL1;
    
    elseif ( strcmp(seq_type, 'RESL') || strcmp(seq_type, 'CSL') )
        SL1 = SL_get_SL_pulse(system, tSL/2, fSL, phase1, SL.freqOffset);
        SL2 = SL_get_SL_pulse(system, tSL/2, fSL, phase2, SL.freqOffset);
                             
    elseif ( strcmp(seq_type, 'PSCSL') || strcmp(seq_type, 'TBSL') )
        SL1 = SL_get_SL_pulse(system, tSL/4, fSL, phase1, SL.freqOffset);
        SL2 = SL_get_SL_pulse(system, tSL/4, fSL, phase2, SL.freqOffset);
    
    elseif strcmp(seq_type, 'BSL')
        SL1 = SL_get_SL_pulse(system, tSL/4, fSL, phase1, SL.freqOffset);
        SL2 = SL_get_SL_pulse(system, tSL/2, fSL, phase2, SL.freqOffset);
        SL3 = SL1;
    end
    end
    
    % T2 case or tSL<100us or fSL<1e-6: use delays
    if ( strcmp(SL.relax_type(loop_SL),'T2') || tSL<100e-6 || fSL<1e-6 )
        SL1 = [];
        SL2 = [];
        SL3 = [];
        if strcmp(seq_type, 'SSL')
            SL1 = mr.makeDelay(tSL);
        elseif ( strcmp(seq_type, 'RESL') || strcmp(seq_type, 'CSL') || strcmp(seq_type, 'hSSL') )
            SL1 = mr.makeDelay(tSL/2);
            SL2 = mr.makeDelay(tSL/2);
        elseif ( strcmp(seq_type, 'PSCSL') || strcmp(seq_type, 'TBSL') || strcmp(seq_type, 'qSSL') )
            SL1 = mr.makeDelay(tSL/4);
            SL2 = mr.makeDelay(tSL/4);
        elseif strcmp(seq_type, 'BSL')
            SL1 = mr.makeDelay(tSL/4);
            SL2 = mr.makeDelay(tSL/2);
            SL3 = mr.makeDelay(tSL/4);
        end
    end

end

%% ----------------------------------------------------------------------------------------
function duration = SL_get_duration(SL)

    %% temporary SL objects
    exc1     = SL.EXC1;
    exc2     = SL.EXC2;
    sl1      = SL.SL1;
    sl2      = SL.SL2;
    sl3      = SL.SL3;
    rfc1     = SL.RFC1;
    rfc2     = SL.RFC2;
    d1       = SL.d1;
    duration = 0;
    
    % add excitation pulse: tip down
    duration = duration + mr.calcDuration(exc1);
    
    % add spin-lock cluster
    duration = duration + mr.calcDuration(d1);
    
    if strcmp(SL.seq_type, 'SSL')
        duration = duration + mr.calcDuration(sl1);
        duration = duration + mr.calcDuration(d1);
    
    elseif strcmp(SL.seq_type, 'hSSL')
        duration = duration + mr.calcDuration(sl1);
        duration = duration + mr.calcDuration(d1);
        duration = duration + mr.calcDuration(sl2);
        duration = duration + mr.calcDuration(d1);
    
    elseif strcmp(SL.seq_type, 'RESL')
        duration = duration + mr.calcDuration(sl1);
        duration = duration + mr.calcDuration(d1);
        duration = duration + mr.calcDuration(sl2);
        duration = duration + mr.calcDuration(d1);
    
    elseif strcmp(SL.seq_type, 'CSL')
        duration = duration + mr.calcDuration(sl1);
        duration = duration + mr.calcDuration(d1);
        duration = duration + mr.calcDuration(rfc1);
        duration = duration + mr.calcDuration(d1);
        duration = duration + mr.calcDuration(sl2);
        duration = duration + mr.calcDuration(d1);
    
    elseif strcmp(SL.seq_type, 'PSCSL')
        duration = duration + mr.calcDuration(sl1);
        duration = duration + mr.calcDuration(d1);
        duration = duration + mr.calcDuration(sl2);
        duration = duration + mr.calcDuration(d1);   
        duration = duration + mr.calcDuration(rfc1);
        duration = duration + mr.calcDuration(d1);
        duration = duration + mr.calcDuration(sl1);
        duration = duration + mr.calcDuration(d1);   
        duration = duration + mr.calcDuration(sl2);
        duration = duration + mr.calcDuration(d1);  
    
    elseif strcmp(SL.seq_type, 'TBSL')
        duration = duration + mr.calcDuration(sl1);
        duration = duration + mr.calcDuration(d1);
        duration = duration + mr.calcDuration(rfc1);
        duration = duration + mr.calcDuration(d1);
        duration = duration + mr.calcDuration(sl2);
        duration = duration + mr.calcDuration(d1);
        duration = duration + mr.calcDuration(sl1);
        duration = duration + mr.calcDuration(d1);
        duration = duration + mr.calcDuration(rfc2);
        duration = duration + mr.calcDuration(d1);
        duration = duration + mr.calcDuration(sl2);
        duration = duration + mr.calcDuration(d1); 
    
    elseif strcmp(SL.seq_type, 'BSL')
        duration = duration + mr.calcDuration(sl1);
        duration = duration + mr.calcDuration(d1);
        duration = duration + mr.calcDuration(rfc1);
        duration = duration + mr.calcDuration(d1);
        duration = duration + mr.calcDuration(sl2);
        duration = duration + mr.calcDuration(d1);
        duration = duration + mr.calcDuration(rfc2);
        duration = duration + mr.calcDuration(d1);
        duration = duration + mr.calcDuration(sl3);
        duration = duration + mr.calcDuration(d1);
    
    end
    
    % add excitation pulse: tip up
    if ~isempty(exc2)
        duration = duration + mr.calcDuration(exc2);
    end

end

%% ----------------------------------------------------------------------------------------
function rf = SL_get_SL_pulse(system, tSL, fSL, phaseOffset, freqOffset)
    rf = mr.makeBlockPulse( pi, system, ...
                           'Duration', tSL, ...
                           'freqOffset', freqOffset, ...
                           'PhaseOffset', phaseOffset, ...
                           'use', 'preparation' );    
    rf.signal = rf.signal*0 + fSL;
end

%% ----------------------------------------------------------------------------------------
function rf = SL_get_comp_pulse(system, tau, phaseOffset)

    % calc amplitude for refocusing pulse
    f1 = 1 / tau;
    
    % create block waveform for composite refocusing pulse
    N          = round(tau/system.rfRasterTime/4)*4;
    comp_amp   = ones(N,1) * f1;
    comp_phase = ones(N,1) * (-pi/2);
    comp_phase(N/4+1 : 3*N/4) = 0; % (-y +x -y) -> (+x)
    signal     = comp_amp .* exp(1i * comp_phase);
    
    % create pulseq rf struct from dummy
    rf             = mr.makeSincPulse( pi, system, 'Duration', N*system.rfRasterTime, 'timeBwProduct', 2, 'use', 'refocusing' );
    rf.signal      = signal;
    rf.phaseOffset = phaseOffset;

end

%% ----------------------------------------------------------------------------------------
function rf = SL_get_block_pulse(system, tau, fa, phaseOffset, freqOffset, use)
rf = mr.makeBlockPulse( fa, system, ...
                       'Duration', tau, ...
                       'freqOffset', freqOffset, ...
                       'PhaseOffset', phaseOffset, ...
                       'use', use );
end