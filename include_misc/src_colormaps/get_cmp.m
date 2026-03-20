function [cmp, xClip] = get_cmp(cmp_name, N, flag_colorLogRemap, cmp_lims, x)
    
    % ----- versions: -----
    % V1: 06.08.2025; M. Gram; University of Wuerzburg

    % ----- input -----
    % cmp_name:           e.g. 'lipari' 'navia' -> source: zenodo.org/records/8409685d + some additionals :D
    % N:                  [ ] size of colormap; use iterpolation of original cmps files
    % flag_colorLogRemap: [0/1] log-like modification according to 10.1002/mrm.30290
    % cmp_lims:           [loLev, upLev] limits of colormap
    % x:                  [Ny x Nx] image to be visualized according to 10.1002/mrm.30290

    % ----- output -----
    % cmp:   [N x 3] colormap
    % xClip: [Ny x Nx] image to be visualized

    % default params
    if nargin<2
        N = [];
    end
    if nargin<3
        flag_colorLogRemap = [];
    end
    if nargin<4
        cmp_lims = [];
    end
    if nargin<5
        x = [];
    end
    if isempty(flag_colorLogRemap)
        flag_colorLogRemap = false;
    end
    if isempty(cmp_lims)
        cmp_lims = [0 1];
    end
    switch cmp_name
        case {'T1', 'R1'}
            cmp_name = 'lipari';
        case {'T2', 'T2*', 'R2', 'R2*', 'T1rho', 'T1ρ', 'T1p', 'R1rho', 'R1ρ', 'R1p'}
            cmp_name = 'navia';
    end

    % load colormap from .mat file
    cmp_files = which('cmp_files.mat');
    try
        temp1 = load(cmp_files, cmp_name);
        temp2 = fieldnames(temp1);
        temp2 = temp2{1};
        eval(['cmp = temp1.' temp2 ';'])
        clear temp1 temp2;
    catch
        error(['can not find colormap: ' cmp_name]);
    end
    
    % interpolate to N
    if ~isempty(N)
        if N ~= size(cmp,1)
                n      = size(cmp,1);
                temp_r = cmp(:,1);
                temp_g = cmp(:,2);
                temp_b = cmp(:,3);
                temp_r = interp1(linspace(0,1,n), temp_r, linspace(0,1,N))';
                temp_g = interp1(linspace(0,1,n), temp_g, linspace(0,1,N))';
                temp_b = interp1(linspace(0,1,n), temp_b, linspace(0,1,N))';
                cmp  = [temp_r, temp_g, temp_b];
                clear n temp_r tempg_g temp_b;
        end
    end
    
    % Fuderer M, et al.
    % Color-map recommendation for MR relaxometry maps.
    % Magn Reson Med. 2025 Feb;93(2):490-506.
    % doi: 10.1002/mrm.30290.
    xClip = [];
    if flag_colorLogRemap==1
        % lookup of the original color map table according to a "log-like" curve
        loLev = cmp_lims(1);
        upLev = cmp_lims(2);
        cmp   = colorLogRemap(cmp, loLev, upLev);

        % modification of the image to be displayed
        if ~isempty(x)   
            eps = (upLev - loLev) / size(cmp, 1);                    
            if loLev < 0
                xClip = arrayfun(@(p) (p < eps) * (loLev - eps) + (p >= eps) * p, x);
            else
                xClip = arrayfun(@(p) (p < eps) * (p < loLev + eps) * (loLev - eps) + (p < eps) * (p >= loLev + eps) * (loLev + eps) + (p >= eps) * p, x);
            end
        end
    end
    
end


%% ------------------------------------------------------------------
% source: https://zenodo.org/records/13142174
% version: 06.08.2025

% colorLogRemap: lookup of the original color map table according to a "log-like" curve.
%   The log-like curve contains a linear part and a logarithmic part; the size of the parts
%   depends on the range (loLev,upLev) 
%
%   Arguments:
%       oriCmap     original colormap, provided as a N*3 matrix
%       loLev       lower level of the range to be displayed
%       upLev       upper level of the range to be displayed
%   Returns:  modified colormap
function logCmap = colorLogRemap(oriCmap, loLev, upLev)
    assert(upLev > 0, 'upper level must be positive');
    assert(upLev > loLev, 'upper level must be larger than lower level');
    
    mapLength = size(oriCmap, 1);
    eInv = exp(-1.0);
    aVal = eInv * upLev;
    mVal = max(aVal, loLev);
    bVal = (1.0 / mapLength) + (aVal >= loLev) * ((aVal - loLev) / (2 * aVal - loLev));
    bVal = bVal+0.0000001;   % This is to ensure that after some math, we get a figure that rounds to 1 ("darkest valid color")
                        % rather than to 0 (invalid color). Note that bVal has no units, so 1E-7 is always a small number    
    logCmap = zeros(size(oriCmap));
    logCmap(1, :) = oriCmap(1, :);
    
    logPortion = 1.0 / (log(mVal) - log(upLev));

    for g = 2:mapLength
        f = 0.0;
        x = g * (upLev - loLev) / mapLength + loLev;
        
        if x > mVal
            % logarithmic segment of the curve
            f = mapLength * ((log(mVal) - log(x)) * logPortion * (1 - bVal) + bVal);
        else
            if (loLev < aVal) && (x > loLev)
                % linear segment of the curve
                f = mapLength * ((x - loLev) / (aVal - loLev) * (bVal - (1.0 / mapLength))) + 1.0;
            end
            
            if (x <= loLev) 
                % lowest valid color
                f = 1.0;
            end
        end
        
        % lookup from original color map
        logCmap(g, :) = oriCmap(min(mapLength, 1 + floor(f)), :);
    end
    
    % Return modified colormap
end
