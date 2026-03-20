clear
pulseq_init();
clearvars -except system

% input parameters example 1:
SPI.geo.design        = 'spiral';
SPI.geo.design_fun    = 'hargreaves';
SPI.geo.sym_rewinder  = false;
SPI.geo.N_interleaves = 24;
SPI.geo.FOV_coeff     = [1 -0.5] * 0.256;
SPI.geo.kmax          = 256 / (2*0.256);
SPI.geo.Ns            = 1e5;
SPI.geo.ds            = 1e-3;
SPI.geo.flag_rv       = 0;

% input parameters example 2:
% SPI.geo.design        = 'spiral';
% SPI.geo.design_fun    = 'log';
% SPI.geo.sym_rewinder  = false;
% SPI.geo.N_loops       = 5;
% SPI.geo.log_coeff     = 1e-6;
% SPI.geo.kmax          = 256 / (2*0.256);
% SPI.geo.Ns            = 1e5;
% SPI.geo.ds            = 1e-3;
% SPI.geo.flag_rv       = 0;

% limit read gradients
SPI.geo.lim_grad = 1.0;
SPI.geo.lim_slew = 1.0;

% calc geometry object
SPI.geo = SPI_minTimeSpiral(SPI.geo, system, 1);
