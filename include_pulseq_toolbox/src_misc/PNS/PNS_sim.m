function [ok, pns_norm, pns_comp]=PNS_sim(obj, slice_orientation)
% calculate PNS using safe model implementation by Szczepankiewicz and Witzel
% assumes safe_pns_prediction package has been downloaded and installed in 
% Matlab path. See http://github.com/filip-szczepankiewicz/safe_pns_prediction
% returns pns levels due to respective axes (normalized to 1 and not to 100%)
% hardware specifications: see safe_example_hw() from
%                          the safe_pns_prediction package. Alternatively a text file
%                          in the .asc format (Siemens) can be passed, e.g. for Prisma
%                          it is MP_GPA_K2309_2250V_951A_AS82.asc (we leave it as an
%                          exercise to the interested user to find were these files
%                          can be acquired from);

% Maximilian Gram V1: 27.09.2024: copied the function seq.calcPNS() from the
% original Pulseq package. extended the option to rotate the slice. removed
% the option to disable plots. set cardiac stimulation to false. the asc
% file for the hardware specification is loaded via the system struct of the
% seq object

% inputs:
%  obj:               the pulseq seq object
%  slice_orientation: 'axial', 'coronal' 'sagittal'
%                     or 3x1 vector containing angles for x, y, z slice rotations

doPlots  = true;
calcCNS  = false;
hardware = obj.sys.ascfile;

% acquire the entire gradient wave form
gw=obj.waveforms_and_times();

% find beginning and end times and resample GWs to a regular sampling raster
tf=[];
tl=[];
for i=1:3
    if size(gw{i},2)>0
        tf(end+1)=gw{i}(1,1);
        tl(end+1)=gw{i}(1,end);
    end
end
nt_min=floor(min(tf)/obj.gradRasterTime+eps); 
nt_max=ceil(max(tl)/obj.gradRasterTime-eps);

% shift raster positions to the centers of the raster periods
nt_min = nt_min + 0.5;
nt_max = nt_max - 0.5;
if (nt_min<0.5)
    nt_min=0.5
end
t_axis=(nt_min:nt_max)*obj.gradRasterTime;
gwr=zeros(length(t_axis),3);
for i=1:3
    if size(gw{i},2)>0
        gwr(:,i)=interp1(gw{i}(1,:),gw{i}(2,:),t_axis,'linear',0);
    end
end

% M. Gram: calculate 3x3 matrix for slice rotation
if nargin<2
    slice_orientation = 'axial';
end
if isstr(slice_orientation)
    if strcmp(slice_orientation, 'axial')
        gwr = ([1 0 0; 0 1 0; 0 0 1] * gwr')';
    elseif strcmp(slice_orientation, 'coronal')
        gwr = ([1 0 0; 0 0 1; 0 1 0] * gwr')';
    elseif strcmp(slice_orientation, 'sagittal')
        gwr = ([0 0 1; 0 1 0; 1 0 0] * gwr')';
    elseif strcmp(slice_orientation, 'all')
        gwr = [ ([1 0 0; 0 1 0; 0 0 1] * gwr')'; ...
                zeros(round(1/obj.gradRasterTime), 3); ...
                ([1 0 0; 0 0 1; 0 1 0] * gwr')'; ...
                zeros(round(1/obj.gradRasterTime), 3); ...
                ([0 0 1; 0 1 0; 1 0 0] * gwr')' ];
        t_axis = (1:size(gwr,1)) * obj.gradRasterTime;
    else
        error('unknown slice orientation');
    end
else
    a = slice_orientation(1);
    b = slice_orientation(2);
    c = slice_orientation(3);
    rot_mat = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)] * [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)] * [cos(c) -sin(c) 0; sin(c) cos(c) 0; 0 0 1];
    gwr     = (rot_mat * gwr')';
    clear a b c;
    slice_orientation = 'oblique';
end

if ischar(hardware)
    % this loads the parameters from the provided text file
    asc=mr.Siemens.readasc(hardware);
    hardware=asc_to_hw(asc,calcCNS);
end

% use the Szczepankiewicz' and Witzel's implementation
[pns_comp,res]=safe_gwf_to_pns(gwr/obj.sys.gamma, NaN*ones(length(t_axis),1), obj.gradRasterTime, hardware); % the RF vector is unused in the code inside but it is zeropaded and exported ... 
% use the exported RF vector to detect and undo zerpopadding
pns_comp=0.01*pns_comp(~isfinite(res.rf),:)';
% calc pns_norm and the final ok/not_ok
pns_norm=vecnorm(pns_comp);
ok=all(pns_norm<1);
% ready
if doPlots
    % plot results
    safe_plot(pns_comp'*100, gwr, obj.gradRasterTime, slice_orientation);
end

end

% local utility functions

function hw = asc_to_hw(asc,useCNS)
% function hw = asc_to_hw(asc)
%
% SAFE model parameters for the asc structure as read from the asc file.
% See comments for units.
% 
% Maxim Zaitsev 08/10/2019

if isfield(asc,'asCOMP') && isfield(asc.asCOMP,'tName')
    hw.name          = asc.asCOMP(1).tName;
else
    hw.name          = 'unknown';
end
%hw.look_ahead    =  1.0; % MZ: this is not a real hardware parameter but a coefficient, with which the final result is multiplied


if isfield(asc,'flGSWDTauX') % older format .asc file
    pns_struct=asc;
    if useCNS 
        error('provided .asc file does not support cardiac stimulation prediction');
    end
elseif isfield(asc,'GradPatSup') % newer format .asc file (e.g. xa61)
    if useCNS 
        if isfield(asc.GradPatSup.Phys, 'CarNS')
            pns_struct=asc.GradPatSup.Phys.CarNS;
        else
            error('provided .asc file does not support cardiac stimulation prediction');
        end
    else
        pns_struct=asc.GradPatSup.Phys.PNS;
    end
else
    error('unknown .asc file format');
end

hw.x.tau1        = pns_struct.flGSWDTauX(1);  % ms
hw.x.tau2        = pns_struct.flGSWDTauX(2);  % ms
hw.x.tau3        = pns_struct.flGSWDTauX(3);  % ms
hw.x.a1          = pns_struct.flGSWDAX(1);
hw.x.a2          = pns_struct.flGSWDAX(2);
hw.x.a3          = pns_struct.flGSWDAX(3);
hw.x.stim_limit  = pns_struct.flGSWDStimulationLimitX;  % T/m/s
hw.x.stim_thresh = pns_struct.flGSWDStimulationThresholdX;  % T/m/s

hw.y.tau1        = pns_struct.flGSWDTauY(1);  % ms
hw.y.tau2        = pns_struct.flGSWDTauY(2);  % ms
hw.y.tau3        = pns_struct.flGSWDTauY(3);  % ms
hw.y.a1          = pns_struct.flGSWDAY(1);
hw.y.a2          = pns_struct.flGSWDAY(2);
hw.y.a3          = pns_struct.flGSWDAY(3);
hw.y.stim_limit  = pns_struct.flGSWDStimulationLimitY;  % T/m/s
hw.y.stim_thresh = pns_struct.flGSWDStimulationThresholdY;  % T/m/s

hw.z.tau1        = pns_struct.flGSWDTauZ(1);  % ms
hw.z.tau2        = pns_struct.flGSWDTauZ(2);  % ms
hw.z.tau3        = pns_struct.flGSWDTauZ(3);  % ms
hw.z.a1          = pns_struct.flGSWDAZ(1);
hw.z.a2          = pns_struct.flGSWDAZ(2);
hw.z.a3          = pns_struct.flGSWDAZ(3);
hw.z.stim_limit  = pns_struct.flGSWDStimulationLimitZ;  % T/m/s
hw.z.stim_thresh = pns_struct.flGSWDStimulationThresholdZ;  % T/m/s

if isfield (asc, 'asGPAParameters')
    hw.x.g_scale = asc.asGPAParameters(1).sGCParameters.flGScaleFactorX;
    hw.y.g_scale = asc.asGPAParameters(1).sGCParameters.flGScaleFactorY;
    hw.z.g_scale = asc.asGPAParameters(1).sGCParameters.flGScaleFactorZ;
else
    hw.x.g_scale = asc.flGCGScaleFactorX; % assume older ASC files like for Trio-Tim
    hw.y.g_scale = asc.flGCGScaleFactorY;
    hw.z.g_scale = asc.flGCGScaleFactorZ;
end

end


