function geo = SPI_importTraj(geo, system, flag_plot)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%   ------------------
%   INPUTS
%   ------
%   geo:  struct with geoemtry parameters; required fields
%         .kmax  - Maximum k-space radius [1/m]
%         .path  - path to a .mat file which has
%                  option 1) Nx1 complex valued k-space  trajectory "kxy" [1/m]
%                  option 2) Nx1 complex valued gradient trajectory "gxy" [Hz/m]
%                  AND the total duration of the waveform "dur" [s]                          
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

geo.lim_grad = [];
geo.lim_slew = [];

%% load complex gxy(t) or kxy(t) and duration
load(geo.path);
if exist('kxy', 'var')
    kxy = kxy(:); % force Nx1
    gxy = [0; diff(kxy)];
    gxy = gxy / max(abs(gxy));
    if ~exist('dur', 'var')
        error('import of k-space trajectory needs a complex valued "kxy" [1/m] and "dur" [s]');
    end
elseif exist('gxy', 'var')
    gxy = gxy(:); % force Nx1
    gxy = gxy / max(abs(gxy));
    if ~exist('dur', 'var')
        error('import of gradient trajectory needs a complex valued "gxy" [Hz/m] and "dur" [s]');
    end
else
    error('import of k-space or gradient trajectory needs a complex valued "kxy" [1/m] or "gxy" [Hz/m] and "dur" [s]');
end

%% interpolate imported trajectory to grad raster time
gxy(1)   = 0;
gxy(end) = 0;
t        = linspace(0, dur, numel(gxy))';
tq       = (0 : 1 : floor(dur/system.gradRasterTime))' *  system.gradRasterTime;
gxy      = interp1(t, gxy, tq);
gxy(1)   = 0;
gxy(end) = 0;
clear t tq;

%% calculate waveform for pulseq
kxy = cumsum(gxy);
kxy = kxy / max(abs(kxy));
kxy = kxy * geo.kmax;
kxy = [real(kxy), imag(kxy)].';
gxy = mr.traj2grad(kxy, 'RasterTime', system.gradRasterTime);
gxy = gxy.';
g   = [gxy(:,1), gxy(:,2), zeros(size(gxy(:,1)))];
clear gxy kxy;

%% calculate effective k-space trajectory and slew rate
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

load(geo.path); % re-load original waveforms

figure()
subplot(1,3,1)
hold on
if exist('kxy', 'var')
    plot(real(kxy), imag(kxy), 'g-', 'LineWidth', 3.0)
end
plot(k(:,1), k(:,2), 'b-', 'LineWidth', 2.0);
xlabel('kx [1/m]');
ylabel('ky [1/m]');
axis image;
set(gca, 'FontName', 'arial', 'FontWeight', 'bold', 'FontSize', 12);

subplot(1,3,2)
hold on;
viscircles([0,0], system.maxGrad, 'Color', 'r', 'LineWidth', 3.0);
if exist('gxy', 'var')
    plot(real(gxy), imag(gxy), 'g-', 'LineWidth', 3.0)
end
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