% simulate_sensitivities.m
%
% Simulates, thanks to the Biot and Savart law, the coils sensitivities in
% a parallel imaging experiment. The coils are equally ditributed around
% the sample. To avoid redundant computations, some symmetries are taken into
% account.
%
% INPUTS:	structure containing: number of coils Nb_coils, radius of coils R, distance of coils
%               D, FOV size FS and resolution Rs (in pixels)
%
% OUTPUT:  	the same structure filled with the FOV maps of coil sensitivities
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group - EPF Lausanne, 2008-02-28

function coil = simulate_sensitivities(coil)

%% Default parameters

if nargin==0
    coil.null = [];
end

if ~isfield(coil,'type')
    coil.type = 'biot';
end

if ~isfield(coil,'Nb_coils')
    coil.Nb_coils = 4; % number of coils
end

if ~isfield(coil,'res')
    coil.res = [128 128];
end

if length(coil.res)==1
    coil.res =[coil.res coil.res];
end

switch coil.type
    case 'biot'
        coil.param.null = [];
        if ~isfield(coil.param,'FS')
            coil.param.FS = 0.26; % FOV width
        end
        if length(coil.param.FS)==1
            coil.param.FS =[coil.param.FS coil.param.FS];
        end
        if ~isfield(coil.param,'D')
            coil.param.D = 0.13; % Distance between the center of the sample and the center of the coils
        end
        if ~isfield(coil.param,'R')
            coil.param.R = 0.08; % radius of the coils
        end
    case 'homogeneous'
        coil.param = [];
end

%% Computation

switch coil.type
    case 'homogeneous'
        S = ones(coil.res(1),coil.res(2),coil.Nb_coils); % for homogeneous sensitivities
        m = 1;
    case 'biot'
        D_theta = 2*pi/coil.Nb_coils;
        theta = (0:coil.Nb_coils-1)*D_theta;
        if isfield(coil.param,'rand')
            theta = theta + 2*pi/3/coil.Nb_coils*randn(1,coil.Nb_coils);
        end
        
        S = zeros(coil.res(1),coil.res(2),coil.Nb_coils);
        m = zeros(1,coil.Nb_coils);
        
        Dx = coil.param.FS(1)/coil.res(1);
        Dy = coil.param.FS(2)/coil.res(2);
        x0 = -coil.param.FS(1)/2:Dx:coil.param.FS(1)/2-Dx; % samples coordinates that are centered
        y0 = -coil.param.FS(2)/2:Dy:coil.param.FS(2)/2-Dy;
        
        theta_done = [];
        
        for ind_coil = 1:coil.Nb_coils % for each coil
            disp(['computing coil ' num2str(ind_coil) ' sensitivity map']);
            % looking for symmetries along y
            ind1 = find(mod(theta(ind_coil)+theta_done,2*pi)==pi); % comment to disable this option
            % looking for symmetries along x
            ind2 = find(mod(theta(ind_coil)+theta_done,2*pi)==0); % comment to disable this option
            if coil.res(1)==coil.res(2)
                % looking for symmetries of rotation 90?
                ind3 = find(mod(theta(ind_coil)-theta_done,2*pi)==pi/2); % comment to disable this option
            end
            if (~isempty(ind1))
                disp('using a symmetry');
                S(:,:,ind_coil) = conj(S(:,end:-1:1,ind1));
            elseif (~isempty(ind2))
                disp('using a symmetry');
                S(:,:,ind_coil) = -conj(S(end:-1:1,:,ind2));
            elseif (~isempty(ind3))
                disp('using a symmetry');
                S(:,:,ind_coil) = -1i*S(end:-1:1,end:-1:1,ind3).';
            else
                S(:,:,ind_coil) = biot_savart_map_c(x0,y0,coil.param.R,coil.param.D,theta(ind_coil));
            end
            tmp = S(:,:,ind_coil);
            theta_done = [theta_done theta(ind_coil)];
            m(ind_coil) = max(abs(tmp(:)));
        end
end

%% Normalization
m = max(m(:));

%for ind_coil = 1:coil.Nb_coils % for each coil
%    S(:,:,ind_coil) = S(:,:,ind_coil)/m;
%end

coil.sensitivity = S;