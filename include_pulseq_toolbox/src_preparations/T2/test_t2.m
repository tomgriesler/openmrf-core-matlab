%% test T2 preparation
clear
pulseq_init();
mrf_mode = 1;

FOV.Nxy    = 128;
FOV.fov_xy = 240 *1e-3;
FOV.dz     = 5 *1e-3;
FOV_init();

%% T2 params
T2.exc_mode   = 'adiabatic_BIR4';
T2.rfc_dur    = 5 *1e-3;
T2.prep_times = [40 80 40 80 40 80 40 80] * 1e-3;
T2 = T2_init(T2, FOV, system);

for loop_T2 = 1 : T2.n_prep
    T2_add();
    seq.addBlock(mr.makeDelay(0.1));
end
seq.plot();
