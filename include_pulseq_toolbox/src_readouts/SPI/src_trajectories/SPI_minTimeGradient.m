function [C, time, g, s, k, phi, sta, stb] = SPI_minTimeGradient(C, dt, g0, gfin, gmax, smax, gamma, ds, rv, show)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%   ----- Input: -----
%   C       - k-space curve [1/m]
%   dt      - sampling time interval [s]
%   g0      - initial gradient amplitude [Hz/m]
%   gfin    - gradient value at the end of the trajectory [Hz/m]
%   gmax    - maximum gradient [Hz/m]
%   smax    - maximum slew [Hz/m/s]
%   gamma   - gyromagnetic ratio [Hz/T]
%   ds      - step size for ODE integration [1/m]
%   rv      - rotationally invariant or variant solution [0/1]
%   show    - show plots while optimizing [0/1]
%   
%   ----- Output: -----
%   C       - reparametrized k-space curve sampled at dt [1/m]
%   time    - total time to get to the end [s]
%   g       - gradiet waveform [Hz/m]
%   s       - slew rate [Hz/m/s]
%   k       - exact k-space corresponding to gradient g [1/m]
%   phi     - geometry constraints on the amplitude vs. arclength
%   sta     - solution for the forward ODE
%   stb     - solution for the backward ODE

% convert input parameters
C    = C  *1e-2;          % [1/m]    -> [1/cm]
dt   = dt *1e3;           % [s]      -> [ms]
g0   = g0   /gamma *1e2;  % [Hz/m]   -> [Gauss/cm]
gfin = gfin /gamma *1e2;  % [Hz/m]   -> [Gauss/cm]
gmax = gmax /gamma *1e2;  % [Hz/m]   -> [Gauss/cm]
smax = smax /gamma *1e-1; % [Hz/m/s] -> [Gauss/cm/ms]
ds   = ds *1e-2;          % [1/m]    -> [1/cm]

% reduce nominal gradient limits to avoid rounding errors
gmax = gmax * 0.99;
smax = smax * 0.99;

% fast method for designing time-optimal gradient waveforms for arbitrary k-space trajectories
[Copt, time, g, s, k, phi, sta, stb] = SPI_find_time_optimized_gradients(C, rv, g0, gfin, gmax, smax, dt, ds);
if isempty(Copt)
    disp(' ');
    disp('   ... optimizing gradient waveform');
    inputs_all    = [C(:); rv+0; g0; gfin; gmax; smax; dt; ds];
    opt_grad_hash = [pulseq_get_path('SPI_find_time_optimized_gradients') 'opt_grad_' pulseq_get_wave_hash(inputs_all) '.mat'];
    [C, time, g, s, k, phi, sta, stb] = minTimeGradient(C, rv, g0, gfin, gmax, smax, dt, ds, show);
    save(opt_grad_hash, 'C', 'time', 'g', 's', 'k', '-v7.3');
    disp('   done & stored for next time!');
else
    C = Copt;
end

% convert output parameters
C    = C        *1e2;    % [1/cm]        -> [1/m]
time = time     *1e-3;   % [ms]          -> [s]
g    = g *gamma *1e-2;   % [Gauss/cm]    -> [Hz/m]
s    = s *gamma *1e1;    % [Gauss/cm/ms] -> [Hz/m/s]
k    = k   *1e2;         % [1/cm]        -> [1/m]
phi  = phi *1e5;         % [1/cm/ms]     -> [1/m/s]
sta  = sta *1e5;         % [1/cm/ms]     -> [1/m/s]
stb  = stb *1e5;         % [1/cm/ms]     -> [1/m/s]

end
