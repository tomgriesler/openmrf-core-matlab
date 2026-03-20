clear

pulseq_init();

FOV.Nxy      = 64;
FOV.fov_xy   = 256  *1e-3;
FOV.dz       = 5    *1e-3;
FOV.z_offset = 0    *1e-3;
FOV_init();

EPI.pe_enable         = 1;   % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
EPI.ro_os             = 1;   % oversampling factor (in contrast to the product sequence we don't really need it)
EPI.readoutTime       = 4.2 *1e-4;   % this controls the readout bandwidth
EPI.partFourierFactor = 1;   % partial Fourier factor: 1: full sampling 0: start with ky=0

EPI.exc_time = 3 *1e-3;
EPI.exc_tbw  = 4;
EPI.exc_mode = 'sinc';

[EPI, ktraj_adc, ktraj_full, ktraj_reco] = EPI_init(EPI, FOV, system);

for j = 1:5
    EPI_add()
end