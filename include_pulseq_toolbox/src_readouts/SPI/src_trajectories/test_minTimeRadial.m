clear
pulseq_init();
clearvars -except system

% input parameters
SPI.geo.design   = 'radial';
SPI.geo.t_adc    = 3.2 *1e-3;
SPI.geo.kmax     = 256 / (2*0.256);

% limit read gradients
SPI.geo.lim_grad = 1.0;
SPI.geo.lim_slew = 1.0;

% calc geometry object
SPI.geo = SPI_minTimeRadial(SPI.geo, system, 1);