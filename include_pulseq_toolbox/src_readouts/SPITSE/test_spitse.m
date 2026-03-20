%%
clear
pulseq_init();

%%
FOV.Nxy      = 120;
FOV.fov_xy   = 240 *1e-3;
FOV.dz       = 5 *1e-3;
FOV.z_offset = 0 *1e-3;

%%
SPITSE.scanmode   = 'run';
SPITSE.segmode    = 'fix';
SPITSE.spmode     = 'cont';
SPITSE.initmode   = 'no';
SPITSE.accmode    = 'vd';
SPITSE.fatsat     = 'off';
SPITSE.seqvar_mod = 'none';
SPITSE.T2prep     = 'on';
SPITSE.EncMode    = 'linear';
SPITSE.plotflag   = '111';
SPITSE.dispflag   = 0;

SPITSE.Trec       = 10 *1e-3;
SPITSE.TE         = 10 *1e-3;
SPITSE.NEcho      = 2;
SPITSE.tEX        = 2.5 *1e-3;
SPITSE.tRef       = 2.0 *1e-3;
SPITSE.flipref    = 60;
SPITSE.flipflag   = 2;
SPITSE.rfex_phase = 0;
SPITSE.tbw        = 4;
SPITSE.slewfac    = 0.99;
SPITSE.gradfac    = 0.99;

[SPITSE, ktraj_adc, ktraj_full, ktraj_reco] = SPITSE_init(SPITSE, FOV, system);

%%

[seq] = SPITSE_add(seq, system, FOV, SPITSE, 1);
seq.plot()
