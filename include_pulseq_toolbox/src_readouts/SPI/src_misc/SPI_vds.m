function [k, g, s, time, r, theta] = SPI_vds(smax, gmax, dt, N, Fcoeff, rmax, gamma)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%   Inputs (Pulseq/SI):
%     smax   [Hz/m/s]   max slew rate (gamma * dB/dt)
%     gmax   [Hz/m]     max gradient amplitude (gamma * B)
%     dt     [s]        sampling period for gradient AND acquisition
%     N      [-]        number of interleaves
%     Fcoeff [m]        FOV polynomial coefficients (see vds.m)
%     rmax   [1/m]      k-space radius at which to stop (typically 1/(2*res))
%     gamma  [Hz/T]     gyromagnetic ratio (default: 42.576e6 for 1H)
%
%   Outputs (Pulseq/SI):
%     k      [1/m]      k-space trajectory (kx + 1i*ky)
%     g      [Hz/m]     gradient waveform (Gx + 1i*Gy)
%     s      [Hz/m/s]   slew waveform (d/dt of g)
%     time   [s]        time vector
%     r      [1/m]      radius vs time
%     theta  [rad]      angle vs time
%
%   Notes:
%   - vds.m internally uses gamma = 4258 Hz/G (hardcoded, i.e. 1H).
%     This wrapper keeps vds.m unchanged and only converts units. If you
%     need a different nucleus, vds.m must be adapted.

if nargin < 7 || isempty(gamma)
    gamma = 42.576e6; % [Hz/T] (1H)
end

% Convert SI -> vds units
% vds expects: gmax [G/cm], smax [G/cm/s], rmax [1/cm], Fcoeff [cm]
Fcoeff = Fcoeff * 1e2;          % [m]      -> [cm]
rmax   = rmax   * 1e-2;         % [1/m]    -> [1/cm]
gmax   = gmax   / gamma * 1e2;  % [Hz/m]   -> [G/cm]
smax   = smax   / gamma * 1e2;  % [Hz/m/s] -> [G/cm/s]

% Call original implementation
[k, g, s, time, r, theta] = vds(smax, gmax, dt, N, Fcoeff, rmax);

% Convert vds units -> SI and column vectors
k = k(:) * 1e2;           % [1/cm]    -> [1/m]
r = r(:) * 1e2;           % [1/cm]    -> [1/m]
g = g(:) * gamma * 1e-2;  % [G/cm]    -> [Hz/m]
s = s(:) * gamma * 1e-2;  % [G/cm/s]  -> [Hz/m/s]
time  = time(:);
theta = theta(:);

end

