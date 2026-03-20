function geo = SPI_minTimeRadialHalf(geo, system, flag_plot)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%   This function designs a radial k-space trajectory with only a half projection.
%   ------------------
%   INPUTS
%   ------
%   geo:  struct with radial geoemtry parameters; required fields
%         .kmax          - Maximum k-space radius [1/m]
%         .design_fun    - option 1: 'minTimeTrap'
%                          needs no further inputs, t_adc is minimized
%         .design_fun    - option 2: 'fixedTimeTrap'
%         .t_adc         - ramp time + flat time
%         .design_fun    - option 3: 'import'
%         .k             - 1d array with k(t)
%         .t             - 1d array with t
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
geo.lim_slew = min([geo.lim_slew, 0.98]);

%% calculate basic radial trajectory g(t)
switch geo.design_fun
    case 'minTimeTrap'
        temp_factor = 0.1;
        temp_check  = 0;
        while temp_check==0
            gx_read = mr.makeTrapezoid('x', 'Area', geo.kmax*temp_factor, 'maxGrad', system.maxGrad*geo.lim_grad, 'maxSlew', system.maxSlew*geo.lim_slew, 'system', system);
            temp_a  = gx_read.flatArea + gx_read.riseTime * gx_read.amplitude / 2;
            if temp_a > geo.kmax * 1.01 % let the read gradient do a little overshoot over kmax
                temp_check = 1;
            else
                temp_factor = temp_factor + 0.001;
            end
        end
        gx_read   = gradTo1d(gx_read,  system);
        gx        = [gx_read(1:end-1); -gx_read];
        g         = [gx, gx*0, gx*0];
        k         = cumsum(g)*system.gradRasterTime;
        s         = diff(g)/system.gradRasterTime;
        geo.t_adc = size(g,1) * system.gradRasterTime;
        clear gx_read gx temp_factor temp_check temp_a;
    case 'fixedTimeTrap'
        gx_read = mr.makeTrapezoid('x', 'FlatArea', geo.kmax, 'FlatTime', geo.t_adc, 'maxGrad', system.maxGrad*geo.lim_grad, 'maxSlew', system.maxSlew*geo.lim_slew, 'system', system);
        gx_read = gradTo1d(gx_read,  system);
        gx      = [gx_read(1:end-1); -gx_read];
        g       = [gx, gx*0, gx*0];
        k       = cumsum(g)*system.gradRasterTime;
        s       = diff(g)/system.gradRasterTime;
        clear gx_read gx temp_factor temp_check temp_a;
    case 'import'
        error('to do. the idea is to import a k(t) trajectory to enable more advanced sampling strategies e.g. with density compensation');
end

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

