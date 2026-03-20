clear
pulseq_init();
clearvars -except system

% note: the key idea of rosette is to spread the fat magnetization in
% transverse plane for the k-space center revists. this will effectively
% attenauate the fat signal. the rosette trajectory is time optimized, but
% reducing the slew rate factor can be used to distribute the fat
% magnetization and to shift the attenuation behaviour.

% input parameters example 1:
SPI.geo.design   = 'rosette';
SPI.geo.N_lobes  = 17;
SPI.geo.kmax     = 256 / (2*0.256);
SPI.geo.Ns       = 1e5;
SPI.geo.ds       = 1e-3;
SPI.geo.flag_rv  = 0;

% input parameters example 2:
% SPI.geo.design   = 'rosette';
% SPI.geo.omega1   = 10.0;
% SPI.geo.omega2   = 4.5;
% SPI.geo.kmax     = 256 / (2*0.256);
% SPI.geo.Ns       = 1e5;
% SPI.geo.ds       = 1e-3;
% SPI.geo.flag_rv  = 0;

% limit read gradients
SPI.geo.lim_grad = 1.0;
SPI.geo.lim_slew = 0.95;

% calc geometry object
SPI.geo = SPI_minTimeRosette(SPI.geo, system, 1);