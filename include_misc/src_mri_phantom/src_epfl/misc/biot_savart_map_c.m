% biot_savart_map_c.cpp
%
% Computes the magnetic field generated on a Cartesian grid by a circular coil.
%
% 
% INPUTS:   x0      vector of x-coordinates of the map samples (M its length)
%           y0      vector of y-coordinates of the map samples (N its length)
%           R       radius of the coil
%           D       distance between the center of the sample and the center of the coil
%           theta   coil's angular position in the sample frame of reference
%
% OUTPUT:  Map of magnetic field complex values (real part: x-component; imaginary part: y-component)
%
% This is a MEX-file for MATLAB.  
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group - EPF Lausanne, 2008-02-28
% revised in october 2010