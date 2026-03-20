function geo = SPI_minTimeRadial(geo, system, flag_plot)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%   This function designs a radial k-space trajectory.
%   ------------------
%   INPUTS
%   ------
%   geo:  struct with radial geoemtry parameters; required fields
%         .kmax          - Maximum k-space radius [1/m]
%         .t_adc         - duration of ADC -> duration of flat time
%   ------
%   geo:  struct with radial geoemtry parameters; optional fields (defaults shown in parentheses):
%         .lim_grad      - Gradient limit scaling factor (1.0)
%         .lim_slew      - Slew-rate limit scaling factor (1.0)
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
if ~isfield(geo, 'lim_grad') || isempty(geo.lim_grad)
    geo.lim_grad = 1.0;
end

if ~isfield(geo, 'lim_slew') || isempty(geo.lim_slew)
    geo.lim_slew = 1.0;
end

% truncate to prevent calculation issues
geo.lim_grad = min([geo.lim_grad, 0.98]);
geo.lim_slew = min([geo.lim_slew, 0.85]);

%% calculate basic radial trajectory g(t)
gx_read = mr.makeTrapezoid('x', 'FlatArea', 2*geo.kmax, 'FlatTime', geo.t_adc, 'maxGrad', system.maxGrad*geo.lim_grad, 'maxSlew', system.maxSlew*geo.lim_slew, 'system', system);
gx_pre  = mr.makeTrapezoid('x', 'Area', -gx_read.area/2, 'maxGrad', system.maxGrad*geo.lim_grad, 'maxSlew', system.maxSlew*geo.lim_slew, 'system', system);
gx_rew  = gx_pre;
gx_read = gradTo1d(gx_read, system);
gx_pre  = gradTo1d(gx_pre,  system);
gx_rew  = gradTo1d(gx_rew,  system);
gx      = [gx_pre(1:end-1); gx_read(1:end-1); gx_rew];
g       = [gx, gx*0, gx*0];
k       = cumsum(g)*system.gradRasterTime;
s       = diff(g)/system.gradRasterTime;
clear gx_read gx_pre gx_rew gx;

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
plot(k(:,1), k(:,2), 'b.', 'LineWidth', 2.0);
xlabel('kx [1/m]');
ylabel('ky [1/m]');
axis image;
set(gca, 'FontName', 'arial', 'FontWeight', 'bold', 'FontSize', 12);

subplot(1,3,2)
hold on;
viscircles([0,0], system.maxGrad, 'Color', 'r', 'LineWidth', 3.0);
plot(g(:,1), g(:,2), 'b.', 'LineWidth', 2.0);
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

function g1d = gradTo1d(g, system)
    g1d = [ ...
           linspace(g.first,      g.amplitude,  g.riseTime/system.gradRasterTime)'; ...
           linspace(g.amplitude,  g.amplitude,  g.flatTime/system.gradRasterTime)'; ...
           linspace(g.amplitude,  g.last,       g.fallTime/system.gradRasterTime)'; ...
          ];
end

