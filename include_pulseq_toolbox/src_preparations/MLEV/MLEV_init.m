function MLEV = MLEV_init(MLEV, FOV, system)

    % Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

    % MLEV parameter initialization
    
    if ~isfield(MLEV, 't_inter')
        MLEV.t_inter = 2*system.gradRasterTime;
    end
    if MLEV.t_inter < 2*system.gradRasterTime
        MLEV.t_inter = 2*system.gradRasterTime;
    end    

    if ~isfield(MLEV, 'ref_phase')
        MLEV.ref_phase = 0;
    end

    if ~isfield(MLEV, 'exc_mode')
        MLEV.exc_mode = 'adiabatic_AHP';
    end

    if strcmp(MLEV.exc_mode, 'adiabatic_AHP')
        if ~isfield(MLEV, 'ahp_exc_time')
            MLEV.ahp_exc_time = 3.0 *1e-3;
        end
        if ~isfield(MLEV, 'ahp_wmax')
            MLEV.ahp_wmax = 600 * 2*pi;
        end
    end

    if strcmp(MLEV.exc_mode, 'adiabatic_BIR4')
        if ~isfield(MLEV, 'bir4_beta')
            MLEV.bir4_beta = 10;            
        end
        if ~isfield(MLEV, 'bir4_kappa')
            MLEV.bir4_kappa = atan(10);            
        end
        if ~isfield(MLEV, 'bir4_dw0')
            MLEV.bir4_dw0 = 30000;            
        end
        if ~isfield(MLEV, 'bir4_f1')
            MLEV.bir4_f1 = 640;            
        end
        if ~isfield(MLEV, 'bir4_alpha')
            MLEV.bir4_alpha = pi/2;            
        end
        if ~isfield(MLEV, 'bir4_dphi')
            MLEV.bir4_dphi = 0.175;            
        end
        if ~isfield(MLEV, 'bir4_tau')
            MLEV.bir4_tau = 0.01;            
        end
    end

    if ~isfield(MLEV, 'crush_nTwists_x')
        MLEV.crush_nTwists_x = 4;   % [] number of 2pi twists in x direction
    end
    if ~isfield(MLEV, 'crush_nTwists_y')
        MLEV.crush_nTwists_y = 4;   % [] number of 2pi twists in y direction
    end
    if ~isfield(MLEV, 'crush_nTwists_z')
        MLEV.crush_nTwists_z = 11.3;   % [] number of 2pi twists in z direction
    end
    if ~isfield(MLEV, 'crush_lim_grad')
        MLEV.crush_lim_grad = 1/sqrt(3); % reduce crusher gradient amplitude from nominal limit
    end    
    if ~isfield(MLEV, 'crush_lim_slew')
        MLEV.crush_lim_slew = 1/sqrt(3); % reduce crusher gradient slew rate from nominal limit to avoid stimulation
    end
    
    MLEV.n_prep         = numel(MLEV.n_mlev);
    MLEV.n_composite    = MLEV.n_mlev * 4;
    MLEV.tau_composite  = 1 / MLEV.fSL;
    MLEV.tau_composite  = round(MLEV.tau_composite/system.gradRasterTime/4) * system.gradRasterTime*4;
    MLEV.t_inter        = ceil(MLEV.t_inter/system.gradRasterTime/2)*system.gradRasterTime*2;
    MLEV.t2p_prep_times = MLEV.n_composite * MLEV.tau_composite; 
    MLEV.t2_prep_times  = MLEV.n_composite * MLEV.t_inter;
    MLEV.fSL            = 1 / MLEV.tau_composite;
    
    % get MLEV objects
    MLEV.MLEV_objs.d1                                                           = mr.makeDelay(system.gradRasterTime);
    MLEV.MLEV_objs.t_inter_1                                                    = mr.makeDelay(MLEV.t_inter/2);
    MLEV.MLEV_objs.t_inter_2                                                    = mr.makeDelay(MLEV.t_inter);
    [MLEV.MLEV_objs.rf_90_td, MLEV.MLEV_objs.rf_90_tu]                          = MLEV_get_objs_EXC(MLEV, system);
    MLEV.MLEV_objs.rf_comp_pos                                                  = MLEV_get_objs_COMP(system, MLEV.tau_composite, MLEV.ref_phase + pi/2);
    MLEV.MLEV_objs.rf_comp_neg                                                  = MLEV_get_objs_COMP(system, MLEV.tau_composite, MLEV.ref_phase - pi/2);
    [MLEV.MLEV_objs.gx_crush, MLEV.MLEV_objs.gy_crush, MLEV.MLEV_objs.gz_crush] = CRUSH_x_y_z(MLEV.crush_nTwists_x, MLEV.crush_nTwists_y, MLEV.crush_nTwists_z, FOV.dx, FOV.dy, FOV.dz, MLEV.crush_lim_grad, MLEV.crush_lim_slew, system);
    [MLEV.MLEV_objs.cycling_list]                                               = MLEV_get_cycling_list(max(MLEV.n_mlev));
    
    % calc total prep durations
    for j = 1 : MLEV.n_prep
        MLEV.duration(j) = MLEV_get_duration(MLEV, j);
    end

end

%% ----------------------------------------------------------------------------------------
function [rf_90_td, rf_90_tu] = MLEV_get_objs_EXC(MLEV, system)

    % AHP: tipdown and tipup
    if strcmp(MLEV.exc_mode, 'adiabatic_AHP')
        AHP      = AHP_get_params(MLEV.ahp_exc_time, MLEV.ahp_wmax);
        rf_90_td = AHP_pulse( system, 'tipdown', AHP.tau, AHP.wmax, AHP.ratio, AHP.beta1, AHP.beta2,  MLEV.ref_phase + pi/2 );
        rf_90_tu = AHP_pulse( system, 'tipup',   AHP.tau, AHP.wmax, AHP.ratio, AHP.beta1, AHP.beta2,  MLEV.ref_phase + pi/2 );
    end
    
    % BIR4: tipdown and tipup
    if strcmp(MLEV.exc_mode, 'adiabatic_BIR4')
        rf_90_td = SIGPY_BIR4(MLEV.bir4_beta, MLEV.bir4_kappa, MLEV.bir4_dw0, MLEV.bir4_f1, MLEV.bir4_alpha, MLEV.bir4_dphi, MLEV.bir4_tau, MLEV.ref_phase + pi/2, 'tipdown', system);
        rf_90_tu = SIGPY_BIR4(MLEV.bir4_beta, MLEV.bir4_kappa, MLEV.bir4_dw0, MLEV.bir4_f1, MLEV.bir4_alpha, MLEV.bir4_dphi, MLEV.bir4_tau, MLEV.ref_phase + pi/2, 'tipup',   system);
    end

end

%% ----------------------------------------------------------------------------------------
function rf = MLEV_get_objs_COMP(system, tau, phaseOffset)

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
function cycling_list = MLEV_get_cycling_list(n_mlev)

    cycling_list = zeros(n_mlev, 4);
    
    for j = 1:n_mlev
        if mod(j,4)==1
            cycling_list(j,:) = [1 1 -1 -1];
        end
        if mod(j,4)==2
            cycling_list(j,:) = [-1 1 1 -1];
        end
        if mod(j,4)==3
            cycling_list(j,:) = [-1 -1 1 1];
        end
        if mod(j,4)==0
            cycling_list(j,:) = [1 -1 -1 1];
        end
    end

end

%% ----------------------------------------------------------------------------------------
function duration = MLEV_get_duration(MLEV, loop_MLEV)

    % temporary MLEV objects
    t_inter_1    = MLEV.MLEV_objs.t_inter_1;
    t_inter_2    = MLEV.MLEV_objs.t_inter_2;
    rf_90_td     = MLEV.MLEV_objs.rf_90_td;
    rf_90_tu     = MLEV.MLEV_objs.rf_90_tu;
    rf_comp_pos  = MLEV.MLEV_objs.rf_comp_pos;
    rf_comp_neg  = MLEV.MLEV_objs.rf_comp_neg;
    cycling_list = MLEV.MLEV_objs.cycling_list';
    cycling_list = cycling_list(:);
    
    % calc MLEV durations
    duration = 0;
    duration = duration + mr.calcDuration(rf_90_td);
    duration = duration + mr.calcDuration(t_inter_1);
    
    for temp_loop = 1 : MLEV.n_composite(loop_MLEV)
    
        if cycling_list(temp_loop) == 1
            duration = duration + mr.calcDuration(rf_comp_pos);
        end
        if cycling_list(temp_loop) == -1
            duration = duration + mr.calcDuration(rf_comp_neg);
        end
        if temp_loop < MLEV.n_composite(loop_MLEV)
            duration = duration + mr.calcDuration(t_inter_2);
        end
    
    end
    clear temp_loop;
    
    duration = duration + mr.calcDuration(t_inter_1);
    duration = duration + mr.calcDuration(rf_90_tu);

end