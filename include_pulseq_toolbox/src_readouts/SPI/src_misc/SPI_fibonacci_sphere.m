function [sph, xyz, area] = SPI_fibonacci_sphere(N, flag_plot)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% Input:
%   N         : number of points to distribute on unit sphere
%   flag_plot : 0 -> off; 1 -> on
%               if 0: no calculation of Voronoi cell areas

% Output:
%   sph  : Spherical coordiantes [az colat] (Nx2)
%   xyz  : Cartesian coordinates [x y z]    (Nx3)
%   area : Voronoi cell area in steradians  (Nx1)

if nargin<2 || isempty(flag_plot)
    flag_plot = 0;
end

%% --- Point generation ---
switch N
    case 4      % Tetrahedron
        xyz = [ ...
             1  1  1;
            -1 -1  1;
            -1  1 -1;
             1 -1 -1];
    case 6      % Octahedron
        xyz = [ ...
             1  0  0;
            -1  0  0;
             0  1  0;
             0 -1  0;
             0  0  1;
             0  0 -1];
    case 8      % Cube
        xyz = [ ...
             1  1  1;
             1  1 -1;
             1 -1  1;
             1 -1 -1;
            -1  1  1;
            -1  1 -1;
            -1 -1  1;
            -1 -1 -1];
    case 12     % Icosahedron
        g = (1 + sqrt(5))/2;
        xyz = [ ...
             0  1  g;
             0 -1  g;
             0  1 -g;
             0 -1 -g;
             1  g  0;
            -1  g  0;
             1 -g  0;
            -1 -g  0;
             g  0  1;
            -g  0  1;
             g  0 -1;
            -g  0 -1];
    case 20     % Dodecahedron
        g = (1 + sqrt(5))/2;
        a = 1/g;
        xyz = [ ...
             1  1  1;
             1  1 -1;
             1 -1  1;
             1 -1 -1;
            -1  1  1;
            -1  1 -1;
            -1 -1  1;
            -1 -1 -1;
             0  a  g;
             0  a -g;
             0 -a  g;
             0 -a -g;
             a  g  0;
             a -g  0;
            -a  g  0;
            -a -g  0;
             g  0  a;
            -g  0  a;
             g  0 -a;
            -g  0 -a];
    otherwise   % Fibonacci sphere
        i  = (0:N-1)';
        ga = pi*(3 - sqrt(5));
        z  = 1 - 2*(i + 0.5)/N;
        r  = sqrt(1 - z.^2);
        phi = i*ga;
        xyz = [r.*cos(phi), r.*sin(phi), z];
end

xyz = xyz ./ vecnorm(xyz,2,2);

%% rotate first point to [0;0;1] -> phi=0, theta=0

% --- Global rotation: force first point to be [0 0 1] ---
v = xyz(1,:).';                 % current first point (3x1)
k = [0; 0; 1];                  % target direction

% If v already on +z (or numerically very close), do nothing
if norm(v - k) > 1e-12

    % Rotation axis (cross product) and angle info
    a = cross(v, k);            % axis proportional (3x1)
    s = norm(a);                % sin(angle)
    c = dot(v, k);              % cos(angle)

    if s < 1e-12
        % v is (anti-)parallel to k.
        % If anti-parallel, rotate 180° around any axis orthogonal to v.
        if c < 0
            % pick an orthogonal axis robustly
            if abs(v(1)) < 0.9
                tmp = [1;0;0];
            else
                tmp = [0;1;0];
            end
            a = cross(v, tmp);
            a = a / norm(a);
            % Rodrigues for 180°: R = -I + 2*a*a'
            R = -eye(3) + 2*(a*a.');
        else
            R = eye(3);
        end
    else
        a = a / s;              % unit rotation axis
        ax = a(1); ay = a(2); az = a(3);

        K = [   0  -az   ay;
               az    0  -ax;
              -ay   ax    0];

        % Rodrigues rotation matrix
        R = eye(3) + K*s + (K*K)*(1 - c);
    end

    % Apply rotation to ALL points
    xyz = (R * xyz.').';

    % Re-normalize for numerical safety
    xyz = xyz ./ vecnorm(xyz,2,2);
end

%% Compute spherical coordinates after rotation
sph = [atan2(xyz(:,2),xyz(:,1)), acos(xyz(:,3))];

area = [];
if flag_plot==0
    return;
end

%% --- Spherical Voronoi via convex hull ---
F = convhull(xyz(:,1), xyz(:,2), xyz(:,3));
M = size(F,1);

Vv = zeros(M,3);
for f = 1:M
    a = xyz(F(f,1),:);
    b = xyz(F(f,2),:);
    c = xyz(F(f,3),:);
    n = cross(b-a,c-a);
    if dot(n,(a+b+c)/3) < 0, n = -n; end
    Vv(f,:) = n / norm(n);
end

adj = cell(N,1);
for f = 1:M
    adj{F(f,1)}(end+1) = f;
    adj{F(f,2)}(end+1) = f;
    adj{F(f,3)}(end+1) = f;
end

area     = nan(N,1);
cellPoly = cell(N,1);

for k = 1:N
    Pk = Vv(adj{k},:);
    if size(Pk,1) < 3, continue; end

    s = xyz(k,:);
    if abs(s(3)) < 0.9, h = [0 0 1]; else, h = [0 1 0]; end
    e1 = cross(h,s); e1 = e1/norm(e1);
    e2 = cross(s,e1);

    ang = atan2(Pk*e2(:), Pk*e1(:));
    [~,ord] = sort(ang);
    Pk = Pk(ord,:);

    cellPoly{k} = Pk;
    area(k) = sph_poly_area(Pk);
end

%% --- Plot (3 subplots) ---
mu   = mean(area,'omitnan');
sig  = std(area,'omitnan');
dA   = area - mu;                         % mean-subtracted areas [sr]
ylim_sym = max(abs(dA(~isnan(dA))));
ylim_sym = max(ylim_sym, 1.2*sig);
ylim_sym = ceil(ylim_sym*10)/10;

figure('Color','w');

subplot(2,2,1); hold on;
[Xs,Ys,Zs] = sphere(100);
surf(Xs,Ys,Zs,'EdgeColor','none','FaceAlpha',0.1);
scatter3(xyz(:,1),xyz(:,2),xyz(:,3), 20, sph(:,1), 'filled');
colormap(get_cmp('vikO', 1000));
axis equal off; view(35,20); camva('manual');
title(sprintf('Points (N=%d)',N));

subplot(2,2,2); hold on;
surf(Xs,Ys,Zs,'EdgeColor','none','FaceAlpha',0.1);
for k = 1:N
    Pk = cellPoly{k};
    if isempty(Pk), continue; end
    Pk = [Pk; Pk(1,:)];
    plot3(Pk(:,1),Pk(:,2),Pk(:,3),'k-','LineWidth',0.8);
end
axis equal off; view(35,20); camva('manual');
title('Voronoi cells');

subplot(2,2,[3,4]); hold on;
plot(dA,'k.','MarkerSize',10);
yline(+sig,'k--','LineWidth',1.0);
yline(-sig,'k--','LineWidth',1.0);
xlim([0 N+1]);
ylim([-ylim_sym ylim_sym]);
grid on;
xlabel('Cell index');
ylabel('\Delta area [sr]');
title(sprintf('\\Delta Voronoi areas (mean-subtracted): std=%.4g sr', sig));

end

function A = sph_poly_area(P)
P  = P ./ vecnorm(P,2,2);
P2 = [P; P(1,:)];
C  = cross(P2(1:end-1,:), P2(2:end,:), 2);
D  = dot(P2(1:end-1,:),   P2(2:end,:), 2);
A  = 2 * atan2(norm(sum(C,1)), 1 + sum(D));
A  = max(0, min(4*pi, A));
end
