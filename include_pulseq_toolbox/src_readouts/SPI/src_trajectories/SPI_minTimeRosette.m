
function geo = SPI_minTimeRosette(geo, system, flag_plot)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%   This function designs a 2D rosette k-space trajectory and computes a time-optimal
%   parametrization subject to gradient and slew-rate constraints using Miki Lustig's
%   minTimeGradient. The resulting gradient, slew-rate, and k-space trajectories are
%   returned in SI/Pulseq units.
%   ------------------
%   NOTE: if flag_plot==1, the attenuation of fat signal (3.5ppm) is
%   calculated and visualized. Use geo.slew_lim to shift the fat signal to
%   a local minimum. Start with geo.slew_lim=1 and reduce ...
%   ------------------
%   INPUTS
%   ------
%   geo:  struct with rosette geoemtry parameters; required fields
%         .kmax          - Maximum k-space radius [1/m]
%         .omega1        - Radial oscillation frequency of the rosette []
%         .omega2        - Angular rotation frequency of the rosette []
%         .N_lobes       - Number of rosette lobes (used if omega1/omega2 empty)
%          If omega1 and omega2 are not specified, they are automatically
%          derived from N_lobes using a golden-ratio-based rule to avoid
%          reducible rosette symmetries.
%   ------
%   geo:  struct with spiral geoemtry parameters; optional fields (defaults shown in parentheses):
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

if ~isfield(geo, 'omega1')
    geo.omega1 = [];
end

if ~isfield(geo, 'omega2')
    geo.omega2 = [];
end

if ~isfield(geo, 'N_lobes')
    geo.N_lobes = [];
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

%% calculate basic rosette trajectory C(p)

if isempty(geo.omega1) || isempty(geo.omega2)
    if isempty(geo.N_lobes)
        error('specify geo.N_lobes or directly set geo.omega1 and geo.omega2');
    end
    if mod(geo.N_lobes, 2)==0
        error('geo.N_lobes should be odd!');
    end
    golden_ratio = (1+sqrt(5))/2;
    geo.omega1   = geo.N_lobes / 2;
    geo.omega2   = round(geo.omega1 / golden_ratio + 0.5) - 0.5;
    geo.addon    = (geo.N_lobes+1) / geo.N_lobes;
else
    if ~isempty(geo.N_lobes)
        warning('geo.N_lobes is ignored if geo.omega1 and geo.omega2 are already specified');
    end
    geo.addon = 1;
end

% design kx + 1i*ky
t   = linspace(0, 2*pi*geo.addon, geo.Ns)';
kxy = geo.kmax * sin(geo.omega1 * t) .* exp(1i * geo.omega2 * t);
kxy = kxy * exp(-1i* angle(kxy(end))); % let the rosette end a (kx,ky)=(kmax,0)

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

%% vis fat attenuation
kx = k(:,1);
ky = k(:,2);

kdist = sqrt(kx.^2+ky.^2);
t0    = find(islocalmin(kdist)) * system.gradRasterTime;
df0   = system.gamma * system.B0 * 3.5 *1e-6;
mxy   = exp(2*pi*1i*t0*df0);
att   = abs(sum(mxy/numel(t0)));
df0_  = 0 : 0.1 : 2.5*df0;
mxy_  = exp(2*pi*1i*t0*df0_);
att_  = abs(sum(mxy_/numel(t0)));

figure();
subplot(2,5,[1 2 3]);
hold on;
plot(df0_, att_*100, 'k-', 'LineWidth', 2);
xline(df0, 'r--');
yline(att*100, 'r--');
plot(df0, att*100, 'ro');
ylabel('attenuation [%]');
title(['attenuation: ' num2str(att*100) '%']);
set(gca, 'FontName', 'arial', 'FontWeight', 'bold', 'FontSize', 12);

subplot(2,5,[6 7 8]);
hold on;
plot(df0_, 20*log10(att_), 'k-', 'LineWidth', 2);
xline(df0, 'r--');
yline(20*log10(att), 'r--');
plot(df0, 20*log10(att), 'ro');
xlabel('Frequency offset [Hz]');
ylabel('attenuation [dB]');
title(['attenuation: ' num2str(20*log10(att)) 'dB']);
set(gca, 'FontName', 'arial', 'FontWeight', 'bold', 'FontSize', 12);

subplot(2,5,[4 5 9 10]);
hold on;
axis image;
for k = 1:numel(mxy)
    plot(mxy(k), '.', 'MarkerSize', 30);
    plot([0 real(mxy(k))], [0 imag(mxy(k))], 'k-');
end
plot(0, 0, 'k.', 'MarkerSize', 40);
xlim([-1 1]*1.1);
ylim([-1 1]*1.1);
xlabel('real(Mxy)');
ylabel('imag(Mxy)');
title('spin evolution in transverse plane');
set(gca, 'FontName', 'arial', 'FontWeight', 'bold', 'FontSize', 12);

sgtitle('change geo.lim_slew to increase attenuation!', 'interpreter', 'none');

end

