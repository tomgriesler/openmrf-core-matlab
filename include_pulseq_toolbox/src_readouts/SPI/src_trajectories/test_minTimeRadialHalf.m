clear
pulseq_init();
clearvars -except system

% input parameters: example 1
SPI.geo.design     = 'radialHalf';
SPI.geo.design_fun = 'minTimeTrap';
SPI.geo.kmax       = 256 / (2*0.256);
SPI.geo.lim_grad   = 1.0;
SPI.geo.lim_slew   = 1.0;

% input parameters: example 2
% SPI.geo.design     = 'radialHalf';
% SPI.geo.design_fun = 'fixedTimeTrap';
% SPI.geo.t_adc      = 3.2 *1e-3;
% SPI.geo.kmax       = 256 / (2*0.256);
% SPI.geo.lim_grad   = 1.0;
% SPI.geo.lim_slew   = 0.85;

% limit read gradients
SPI.geo.lim_grad = 1.0;
SPI.geo.lim_slew = 1.0;

% calc geometry object
SPI.geo = SPI_minTimeRadialHalf(SPI.geo, system, 1);