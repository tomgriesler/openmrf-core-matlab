function phantom = mg_numerical_phantom(Nxy, geo, coil_params)

if nargin<3
    coil_params = [];
end

%% define phantom basics 

% geoemtry similar to the Caliber MRI phantom
if isempty(geo)
    n_inner = 5;    % number of circular phantom inserts (inner)
    n_outer = 9;    % number of circular phantom inserts (outer)
    r_phant = 0.9;  % radius of main phantom
    r_outer = 0.65; % radius of circle along which outer phantoms lay
    r_inner = 0.3;  % radius of circle along which inner phantoms lay
    r_insrt = 0.1;  % radius of each circular insert
else
    n_inner = geo.n_inner;
    n_outer = geo.n_outer;
    r_phant = geo.r_phant;
    r_outer = geo.r_outer;
    r_inner = geo.r_inner;
    r_insrt = geo.r_insrt;
end

% total number of circular phantom inserts
n_phant = n_inner + n_outer;

% init basic circular phantom
circ.FOV = [1, 1];
circ.region{1}.type   = 'ellipse';
circ.region{1}.center = [0, 0];
circ.region{1}.angle  = 0;
circ.region{1}.weight = 1;
circ.region{1}.width  = r_phant * [1, 1];
res                   = Nxy * [1, 1];
mask2D                = imfill(RasterizePhantom(circ, res, 1)==1, 'holes');

%% simulate coils

if ~isempty(coil_params)

    % https://bigwww.epfl.ch/algorithms/mriphantom/
    % https://bigwww.epfl.ch/algorithms/mriphantom/MRIPhantomv0-8.zip
    % Guerquin-Kern M, Lejeune L, Pruessmann KP, Unser M.
    % Realistic analytical phantoms for parallel magnetic resonance imaging.
    % IEEE Trans Med Imaging. 2012 Mar;31(3):626-36.
    % doi: 10.1109/TMI.2011.2174158
    
    Ncoils = coil_params.Nb_coils;
    
    % parameter of the model: polynomial degree or bandwith
    if strcmp(coil_params.sens_model, 'polynomial')
        param = 8;
    end
    if strcmp(coil_params.sens_model, 'sinusoidal')
        param = 6;
    end
    coil_params = simulate_sensitivities(coil_params);
    
    % coil simultion
    cmaps    = zeros(Ncoils, Nxy, Nxy);
    s        = cell(1, Ncoils);
    sens_num = cell(1, Ncoils);
    residue  = cell(1, Ncoils);
    for j = 1 : Ncoils
        sens             = coil_params.sensitivity(:,:,j);
        sens             = sens / max(sens(mask2D));
        s{j}             = SensFitting(sens, coil_params.sens_model, param, mask2D);
        [~, sens_num{j}] = RasterizePhantom(circ, res, s{j});
        residue{j}       = sens_num{j}(mask2D) - sens(mask2D);
        cmaps(j,:,:)     = sens_num{j} .* mask2D;
        clear sens;
    end
    
    cmaps = cmaps / max(abs(cmaps(:)));

else
    % load measured cmaps from NIST scan
    load('cmaps_20_330_330_NIST.mat'); % measured @3T CimaX
    Ncoils = size(cmaps, 1);
    cmaps = interp_cmaps_NIST(cmaps, max(sum(mask2D))+2, Nxy);
    cmaps = cmaps .* permute(repmat(mask2D,1,1,Ncoils), [3,1,2]);
end

%% create circular phantom inserts -> ind_map 

% angles of circular inserts
phi_inner = linspace(0, 2*pi, n_inner+1); phi_inner(end) = []; % inner circle
phi_outer = linspace(0, 2*pi, n_outer+1); phi_outer(end) = []; % outer circle
phi_inner = phi_inner + pi/3; % twist inner circle w.r.t. outer circle
phi_outer = phi_outer + pi/5; % twist outer circle w.r.t. inner circle

% inner circles
ind_map  = zeros(Nxy,Nxy);
xC_inner = r_inner * cos(phi_inner); xC_inner = xC_inner*(Nxy/2) + Nxy/2;
yC_inner = r_inner * sin(phi_inner); yC_inner = -yC_inner*(Nxy/2) + Nxy/2;
rC_inner = r_insrt * ones(length(xC_inner),1)'; rC_inner = rC_inner * (Nxy/2);
circles_inner = transpose([xC_inner; yC_inner; rC_inner]);
for j=1:length(circles_inner)
    for y=1:Nxy
    for x=1:Nxy
        if ( (x-circles_inner(j,1))^2 + (y-circles_inner(j,2))^2 ...
                < circles_inner(j,3)^2)
            ind_map(y,x) = j+1;
        end
    end
    end
end
clear j x y

% outer circles 
xC_outer = r_outer * cos(phi_outer); xC_outer =  xC_outer*(Nxy/2) + Nxy/2;
yC_outer = r_outer * sin(phi_outer); yC_outer = -yC_outer*(Nxy/2) + Nxy/2;
rC_outer = r_insrt * ones(length(xC_outer),1)'; rC_outer = rC_outer * (Nxy/2);
circles_outer = transpose([xC_outer; yC_outer; rC_outer]);
for j=1:length(circles_outer)
    for y=1:Nxy
    for x=1:Nxy
        if ( (x-circles_outer(j,1))^2 + (y-circles_outer(j,2))^2 ...
                < circles_outer(j,3)^2)
            ind_map(y,x) = j+n_inner+1;
        end
    end
    end
end
clear j x y

% background -> index 1
ind_map((mask2D==1) & (ind_map==0)) = 1;

% useful for converting 1d <-> 2d
pos_map = reshape(1:Nxy^2, [Nxy, Nxy]);

%% output phantom as struct
phantom.Nxy         = Nxy;
phantom.Ncoils      = Ncoils;
phantom.n_phant     = n_phant;
phantom.geo.n_inner = n_inner;
phantom.geo.n_outer = n_outer;
phantom.geo.r_phant = r_phant;
phantom.geo.r_inner = r_inner;
phantom.geo.r_outer = r_outer;
phantom.geo.r_insrt = r_insrt;
phantom.coil_params = coil_params;
phantom.cmaps       = cmaps;
phantom.ind_map     = ind_map;
phantom.pos_map     = pos_map;
phantom.mask2D      = mask2D;
phantom.mask1D      = mask2D(:)==1;

end

%% interpolation function for measured cmaps
function cmaps_out = interp_cmaps_NIST(cmaps, m, M)
    [Ncoils, n, ~] = size(cmaps);
    x        = linspace(-0.5, 0.5, n);
    xq       = linspace(-0.5, 0.5, m);
    [X, Y]   = ndgrid(x, x);
    [Xq, Yq] = ndgrid(xq, xq);
    
    cmaps_interp = complex(zeros(Ncoils, m, m));
    for c = 1:Ncoils
        Sc = squeeze(cmaps(c,:,:));
        Fr = griddedInterpolant(X,Y,real(Sc), 'linear', 'nearest');
        Fi = griddedInterpolant(X,Y,imag(Sc), 'linear', 'nearest');
        cmaps_interp(c,:,:) = complex(Fr(Xq,Yq), Fi(Xq,Yq));
    end

    start     = floor((M - m)/2) + 1;
    stop      = start + m - 1;
    cmaps_out = zeros(Ncoils, M, M);
    cmaps_out(:, start:stop, start:stop) = cmaps_interp;
end
