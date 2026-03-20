clear
pulseq_init();

%% FOV geometry
FOV.Nx       = 256;
FOV.Ny       = 256;
FOV.fov_x    = 240 *1e-3;
FOV.fov_y    = 240 *1e-3;
FOV.dz       = 5 *1e-3;
FOV.z_offset = 0 *1e-3;
FOV_init();

%% TSE sequence parameters
TSE.n_echo    = 4;
TSE.TE        = 12   *1e-3;
TSE.TR        = 500   *1e-3;
TSE.Ndummy    = 2;
TSE.exc_time  = 2.5 *1e-3;
TSE.rfc_time  = 2.5 *1e-3;
TSE.exc_tbw   = 4;
TSE.rfc_tbw   = 4;
TSE.t_acq     = 6.4 *1e-3 + 2*system.adcDeadTime;
TSE.os_mode   = 1;  % read oversampling: 0 off, 1 on
TSE.mode_exc  = 'sinc'; % 'sigpy_SLR' or 'sinc'
TSE.mode_rfc  = 'sinc'; % 'sigpy_SLR' or 'sinc'
TSE.enc_mode  = 'centric';

[TSE, ktraj_adc, ktraj_full] = TSE_init(TSE, FOV, system);

%%
for loop_TR = 1-TSE.Ndummy : TSE.nex
    TSE_add();
end

seq.plot()
