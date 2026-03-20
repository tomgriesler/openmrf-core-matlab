function geo = SPI_minTimeSpiral(geo, system, flag_plot)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%   This function designs a 2D spiral k-space trajectory and computes a time-optimal
%   parametrization subject to gradient and slew-rate constraints using Miki Lustig's
%   minTimeGradient. The resulting gradient, slew-rate, and k-space trajectories are
%   returned in SI/Pulseq units.
%   ------------------
%   INPUTS
%   ------
%   geo:  struct with spiral geoemtry parameters; required fields
%         .kmax          - Maximum k-space radius [1/m]
%         .design_fun    - option 1: 'hargreaves'
%         .N_interleaves - number of spiral interleveaves for sampling the k-space center
%         .FOV_coeff     - FOV*[1 0]    -> archimediean
%                          FOV*[1 -0.5] -> typical choice for MRF
%                          FOV*[1 -1]   -> logarithmic
%         .design_fun    - option 2: 'log'
%         .N_loops       - number of spiral loops
%         .log_coeff     - geometry coefficient; 0...1 -> archimediean...logarithmic...twirl
%   ------
%   geo:  struct with spiral geoemtry parameters; optional fields (defaults shown in parentheses):
%         .alpha         - Bezier control length along end tangent (default: 0.2*|k_end|)
%         .beta          - Bezier control length along radial direction (default: 0.2*|k_end|)
%         .lim_grad      - Gradient limit scaling factor (1.0)
%         .lim_slew      - Slew-rate limit scaling factor (1.0)
%         .Ns            - Number of samples for initial curve parametrization (1e5)
%         .ds            - Arc-length step for minTimeGradient (1e-3)
%         .flag_rv       - Rotationally invariant solver flag (0)
%   ------
%   system: struct with scanner hardware parameters (Pulseq conventions)
%           .maxGrad [Hz/m]
%           .maxSlew [Hz/m/s]
%           .adcDeadTime [s]
%           .adcRasterTime [s]
%           .gradRasterTime [s]
%           .gamma [Hz/T]
%           .B0 [T]
%   ------
%   flag_plot : (optional) visualization flag
%       0 or omitted  - no plots
%       1             - plot k-space, gradient, and slew rate
%   ------------------
%   OUTPUT
%   -------
%   geo:  final struct with spiral geoemtry parameters and waveforms
%         .g : [Nt x 3] gradient waveform [Hz/m]
%              Columns correspond to gx, gy, gz.
%
%         .s : [Nt x 3] slew-rate waveform [Hz/m/s]
%              Time derivative of g.
%
%         .k : [Nt x 3] k-space trajectory [1/m]
%              Time integral of g.

%% set defaults

if ~isfield(geo, 'alpha')
    geo.alpha = [];
end

if ~isfield(geo, 'beta')
    geo.beta = [];
end

if ~isfield(geo, 'sym_rewinder') || isempty(geo.sym_rewinder)
    geo.sym_rewinder = false;
end

if ~isfield(geo, 'sample_rewinder') || isempty(geo.sample_rewinder)
    geo.sample_rewinder = false;
end

if ~isfield(geo, 'lim_grad') || isempty(geo.lim_grad)
    geo.lim_grad = 1.0;
end

if ~isfield(geo, 'lim_slew') || isempty(geo.lim_slew)
    geo.lim_slew = 1.0;
end

if ~isfield(geo, 'Ns') || isempty(geo.Ns)
    geo.Ns = 1e5;
end

if ~isfield(geo, 'ds') || isempty(geo.ds)
    geo.ds = 1e-3;
end

if ~isfield(geo, 'flag_rv') || isempty(geo.flag_rv)
    geo.flag_rv = 0;
end

%% calculate basic spiral trajectory C(p)

% design kx + 1i*ky
if strcmp(geo.design_fun, 'log') % based on a simple analytical design function
    phi = linspace(0, 2*pi * geo.N_loops, geo.Ns)';
    kr  = geo.kmax * (exp(geo.log_coeff * phi) - 1) / (exp(geo.log_coeff * 2*pi*geo.N_loops) - 1);
    kxy = kr .* (cos(phi) + 1i * sin(phi));
    geo.N_interleaves = [];
    geo.FOV_coeff     = [];
elseif strcmp(geo.design_fun, 'hargreaves') % based on the variable-density toolbox of Brian Hargreaves
    kxy = SPI_vds(system.maxSlew*geo.lim_slew, system.maxGrad*geo.lim_grad, 1e-6, geo.N_interleaves, geo.FOV_coeff, geo.kmax);
    geo.N_loops   = [];
    geo.log_coeff = [];
end
kxy = kxy * exp(-1i* angle(kxy(end))); % let the spiral end a (kx,ky)=(kmax,0)

% design kxy rewinder
if geo.sym_rewinder
    kxy = [kxy; flipud(conj(kxy))];
else
    Nrew = round(numel(kxy) * 0.1);
    kxy  = SPI_addBezierRewinder(kxy, geo.alpha, geo.beta, Nrew);
end

% combine to Nx3
kx = real(kxy);
ky = imag(kxy);
C  = [kx, ky, zeros(size(kx))];
C(end,1:2) = 0;

% remove duplicate neighbours
tol  = 1e-9;
keep = abs(diff(C)) > tol;
keep = [true; sum(keep,2)>0];
C    = C(keep,:);

%% find fastest parametrization p(t) to run through C(p)
[~, ~, g] = SPI_minTimeGradient( C, ...
                                 system.gradRasterTime, ...
                                 0, ...
                                 0, ...
                                 system.maxGrad * geo.lim_grad, ...
                                 system.maxSlew * geo.lim_slew, ...
                                 system.gamma, ...
                                 geo.ds, ...
                                 geo.flag_rv, ...
                                 0);

k = cumsum(g) * system.gradRasterTime;
s = diff(g)   / system.gradRasterTime;

%% export waveforms to geo struct
geo.g = g;
geo.k = k;
geo.s = s;

%% vis trajectory
if nargin<3 || flag_plot==0
    return;
end

figure()
subplot(1,3,1)
hold on
plot(kx, ky, 'r-', 'LineWidth', 3.0)
plot(k(:,1), k(:,2), 'b-', 'LineWidth', 2.0);
xlabel('kx [1/m]');
ylabel('ky [1/m]');
axis image;
set(gca, 'FontName', 'arial', 'FontWeight', 'bold', 'FontSize', 12);

subplot(1,3,2)
hold on;
viscircles([0,0], system.maxGrad, 'Color', 'r', 'LineWidth', 3.0);
plot(g(:,1), g(:,2), 'b-', 'LineWidth', 2.0);
xlabel('gx [Hz/m]');
ylabel('gy [Hz/m]');
axis image;
set(gca, 'FontName', 'arial', 'FontWeight', 'bold', 'FontSize', 12);

subplot(1,3,3)
hold on;
viscircles([0,0], system.maxSlew, 'Color', 'r', 'LineWidth', 3.0);
plot(s(:,1), s(:,2), 'b.', 'LineWidth', 2.0);
xlabel('sx [Hz/m/s]');
ylabel('sy [Hz/m/s]');
axis image;
set(gca, 'FontName', 'arial', 'FontWeight', 'bold', 'FontSize', 12);

end

