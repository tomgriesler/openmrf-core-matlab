%% test spin-lock preparation
clear;
pulseq_init();
mrf_mode = 1;

FOV.Nxy    = 128;
FOV.fov_xy = 240 *1e-3;
FOV.dz     = 5 *1e-3;
FOV_init();

%% params
SL.relax_type = {'T1p'};
SL.seq_type   = {'TBSL'};
SL.tSL        = [8 16 24 32 40 48 56 64 72 80] *1e-3;
SL.fSL        = [300 300 300 300 300 300 300 300 300 300];
SL.inv_mode   = {'off'};

SL.freqOffset = 0.0; % [Hz] frequency offset of SL pulses

% excitation pulses
SL.exc_mode       = 'adiabatic_AHP';  % 'adiabatic_AHP', 'adiabatic_BIR4','sinc', 'sigpy_SLR' or 'bp'
SL.exc_flip_angle = pi/2;         % [rad] flip angle
SL.exc_time       = 3.0 *1e-3;    % [s] excitation time
SL.exc_tbw        = 6;            % [ ] time bandwith product, only for sinc/sigpy_SLR
SL.adia_wmax      = 600 * 2*pi;   % [rad/s] amplitude of adiabatic pulse

% refocusing pulses
SL.rfc_mode       = 'sinc';       % 'sinc', 'sigpy_SLR', 'bp' or 'comp'
SL.rfc_flip_angle = pi;           % [rad] flip angle
SL.rfc_time       = 3.0 *1e-3;    % [s] refocusing time
SL.rfc_tbw        = 4;            % [ ] time bandwith product, only for sinc/sigpy_SLR

SL = SL_init(SL, FOV, system);

for loop_SL = 1:SL.nSL
    SL_add();
    seq.addBlock(mr.makeDelay(0.1));
end
seq.plot()
