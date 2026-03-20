%% test mlev preparation
clear
pulseq_init();
mrf_mode = 1;

FOV.Nxy    = 128;
FOV.fov_xy = 240 *1e-3;
FOV.dz     = 5 *1e-3;
FOV_init();

%% mlev params
MLEV.n_mlev   = [ 1 2 3 4 5 6 7 8 9 10]; % number of MLEV4 preps
MLEV.fSL      = 250;             % [Hz] eff spin-lock field strength
MLEV.t_inter  = system.gradRasterTime; % [s]  inter pulse delay for T2 preparation
MLEV.exc_mode = 'adiabatic_AHP'; % 'adiabatic_BIR4' or 'adiabatic_AHP'
MLEV = MLEV_init(MLEV, FOV, system);

for loop_MLEV = 1 : MLEV.n_prep
    MLEV_add();
    seq.addBlock(mr.makeDelay(0.1));
end
seq.plot();
