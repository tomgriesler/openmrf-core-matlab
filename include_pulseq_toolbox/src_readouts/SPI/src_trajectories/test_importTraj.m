clear
pulseq_init();
clearvars -except system

% inputs
SPI.geo.design   = 'import';
SPI.geo.path     = 'test_import_kxy';
% SPI.geo.path     = 'test_import_gxy';
SPI.geo.kmax     = 256 / (2*0.256);

% limit read gradients
SPI.geo.lim_grad = 1.0;
SPI.geo.lim_slew = 1.0;

% calc geometry object
SPI.geo = SPI_importTraj(SPI.geo, system, 1);