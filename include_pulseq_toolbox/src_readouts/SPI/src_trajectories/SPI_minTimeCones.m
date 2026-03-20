function geo = SPI_minTimeCones(geo, system, flag_plot)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%   This function designs a 3D spiral cones k-space trajectory and computes a time-optimal
%   parametrization subject to gradient and slew-rate constraints using Miki Lustig's
%   minTimeGradient. The resulting gradient, slew-rate, and k-space trajectories are
%   returned in SI/Pulseq units.
%   ------------------
%   INPUTS
%   ------
%   geo:  struct with spiral cones geoemtry parameters; required fields
%         .kmax          - Maximum k-space radius [1/m]
%         .N_loops       - number of spiral loops
%         .theta         - cone angle to +kz [rad]
%         .FOV_coeff     - [1 0]    -> archimediean
%                          [1 -0.5] -> typical choice for MRF
%                          [1 -1]   -> logarithmic
%   ------
%   geo:  struct with spiral cones geoemtry parameters; optional fields (defaults shown in parentheses):
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
%   geo:  final struct with spiral cones geoemtry parameters and waveforms
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

if (geo.FOV_coeff(1) + geo.FOV_coeff(2)) <= 0
    error('Invalid geo.FOV_coeff: require a+b > 0 so FOV stays positive up to the periphery.');
end

%% calculate basic spiral cone trajectory C(p)
k    = linspace(0, geo.kmax, geo.Ns).';
kz   = k * cos(geo.theta);
kxy  = k * sin(geo.theta);
a    = geo.FOV_coeff(1);
b    = geo.FOV_coeff(2);
r    = kxy;
u    = r / max(kxy);
den  = (a + b/2);
phi  = 2*pi*geo.N_loops * (a*u + (b/2)*u.^2) / den;
kx   = kxy .* cos(phi);
ky   = kxy .* sin(phi);
C    = [kx ky kz];  % Nx3
C    = SPI_addBezierRewinder3D(C, geo.alpha, geo.beta, round(0.1*geo.Ns));

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

ns      = 256;
[th,ph] = meshgrid(linspace(0,pi,ns), linspace(0,2*pi,ns));
Xs      = sin(th).*cos(ph);
Ys      = sin(th).*sin(ph);
Zs      = cos(th);
gvio    = g;
svio    = s;
gvio( sqrt(g(:,1).^2+g(:,2).^2+g(:,3).^2)<system.maxGrad, : ) = [];
svio( sqrt(s(:,1).^2+s(:,2).^2+s(:,3).^2)<system.maxSlew, : ) = [];

figure('Color','w');
subplot(1,3,1)
hold on;
surf(Xs*geo.kmax, Ys*geo.kmax, Zs*geo.kmax, 'EdgeColor', 'none', 'FaceAlpha', 0.05);
plot3(k(:,1), k(:,2), k(:,3), 'b-', 'LineWidth', 2.0);
xlabel('kx [1/m]');
ylabel('ky [1/m]');
zlabel('kz [1/m]');
axis image;
set(gca, 'FontName', 'arial', 'FontWeight', 'bold', 'FontSize', 12);
view([-37.5 -30]);

subplot(1,3,2)
hold on;
surf(Xs*system.maxGrad, Ys*system.maxGrad, Zs*system.maxGrad, 'EdgeColor', 'none', 'FaceAlpha', 0.05);
plot3(g(:,1), g(:,2), g(:,3), 'b-', 'LineWidth', 2.0);
plot3(gvio(:,1), gvio(:,2), gvio(:,3), 'r-', 'LineWidth', 2.0);
xlabel('gx [Hz/m]');
ylabel('gy [Hz/m]');
zlabel('gz [Hz/m]');
axis image;
set(gca, 'FontName', 'arial', 'FontWeight', 'bold', 'FontSize', 12);
view([-37.5 -30]);

subplot(1,3,3)
hold on;
surf(Xs*system.maxSlew, Ys*system.maxSlew, Zs*system.maxSlew, 'EdgeColor', 'none', 'FaceAlpha', 0.05);
plot3(s(:,1), s(:,2), s(:,3), 'b.', 'LineWidth', 2.0);
plot3(svio(:,1), svio(:,2), svio(:,3), 'r.', 'LineWidth', 2.0);
xlabel('sx [Hz/m/s]');
ylabel('sy [Hz/m/s]');
zlabel('sz [Hz/m/s]');
axis image;
set(gca, 'FontName', 'arial', 'FontWeight', 'bold', 'FontSize', 12);
view([-37.5 -30]);

end

