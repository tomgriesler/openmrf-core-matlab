function DATA = mg_phantom_sim(Mxy, mask2D, cmaps, ktraj, proj_id, T2s, dw0, dwell, fov)

% forward model for numerical phantom simulation
% Maximilian Gram, University of Würzburg, v1, 13.03.2026

% ----- inputs: -----
% Mxy:     NR x Nvoxel, initial magnetization at beginning of adc
% mask2D:  Nxy x Nxy
% cmaps:   NCoils x Nxy x Nxy
% ktraj:   2 x Nunique x NRead, [1/m]
% proj_id: NR x 1

% ----- optional inputs: -----
% T2s:   [s]  Nxy x Nxy, T2* relaxation time
% df0:   [Hz] Nxy x Nxy, offresonance
% dwell: [s]  dwell time of adc
% fov:   [m]  field of view

% ----- output: -----
% DATA: NCoils x NR x Nread

% check simulation mode
sim_mode = 'nufft'; % default: simple forward model based on nufft, no off-resonance, no T2s
if nargin>6 && ~isempty(T2s) && ~isempty(dw0) && ~isempty(dwell) && ~isempty(fov)
    sim_mode = 'phys'; % forward model based on "physical" signal equation
else
    T2s   = [];
    dw0   = [];
    dwell = [];
    fov   = [];
end

% convert inputs
Mxy     = single(Mxy).';
cmaps   = permute(single(cmaps), [2,3,1]);
ktraj   = single(ktraj);
proj_id = single(proj_id);
T2s     = single(T2s);
dw0     = single(dw0);
dwell   = single(dwell);
fov     = single(fov);
mask2D  = mask2D==1;

% get dimensions
[~, NR]          = size(Mxy);
[Nxy, ~, NCoils] = size(cmaps);
NRead            = size(ktraj, 3);

tic;
fprintf('\n');
fprintf(['numerical phantom simulation based on ' sim_mode ' \n']);

% start simulation
switch sim_mode
    case 'nufft'
        % init array for simulated MRI signal
        DATA = zeros(NR, NRead, NCoils, 'single');
        
        % init NUFFT operator
        kx   = single(squeeze(ktraj(1,:,:)))';
        ky   = single(squeeze(ktraj(2,:,:)))';
        kmax = max(sqrt(kx(:).^2 + ky(:).^2));
        kx   = kx/kmax*Nxy/2;    
        ky   = ky/kmax*Nxy/2;
        FT   = parallel.pool.Constant(@() NUFFT(kx/Nxy+1i*ky/Nxy, kx*0+1, [0 0], Nxy*[1 1]) );
        clear ktraj kx ky kmax fov Ni nufft_args G dcf_all;

        % get mask indices
        idx = reshape(1:(Nxy)^2, Nxy,Nxy);
        idx = single(idx(mask2D));
        
        % simple forward model:
        % calc complex images for each TR -> add coils sens -> apply nufft -> apply sampling mask
        parfor j = 1:NR
            temp_mxy    = Mxy(:,j);
            temp_image  = zeros(Nxy, Nxy);
            temp_image(idx) = temp_mxy;
            temp_image  = repmat(temp_image, 1, 1, NCoils) .* cmaps;
            temp_signal = FT.Value * temp_image;
            temp_signal = squeeze(temp_signal(:,proj_id(j),:));
            DATA(j,:,:) = temp_signal;
        end
        DATA = permute(DATA, [3, 1, 2]); % NCoils x NR x NRead

    case 'phys'
        % init array for simulated MRI signal
        DATA = zeros(NRead, NR, NCoils, 'single');
        
        % convert maps to 1D
        dw0 = dw0(mask2D);
        T2s = T2s(mask2D);
        cmaps1D = zeros(numel(dw0), NCoils);
        for j=1:NCoils
            temp = squeeze(cmaps(:,:,j));
            cmaps1D(:,j) = temp(mask2D);
        end
        clear temp cmaps;

        % define real space of phantom in 1D
        [x, y] = meshgrid( -Nxy/2 : Nxy/2-1, -Nxy/2 : Nxy/2-1 );
        x = x';
        y = y';
        x = x(:) * fov/Nxy; % [m]
        y = y(:) * fov/Nxy; % [m]
        x(mask2D(:)==0) = [];
        y(mask2D(:)==0) = [];
        
        % init k-space
        kx  = 2*pi * squeeze(ktraj(1,proj_id,:)); % [rad/m]
        ky  = 2*pi * squeeze(ktraj(2,proj_id,:)); % [rad/m]
        clear ktraj;

        % init off-resonance and decay during sampling
        phi_off   = dwell * repmat(dw0, 1, NR); % [rad]
        t2s_decay = dwell ./ repmat(T2s, 1, NR); % [ ]
        
        % case: NR=1
        if NR==1
            kx = kx.';
            ky = ky.';
        end
               
        % physical forward model:
        % magnetization state before adc -> T2s decay -> df0 dephasing -> gradient dephasing -> apply coil sens
        parfor j=1:NRead
                temp_signal = Mxy .* ...                      % simulated magnetization
                              exp(-t2s_decay*j) .* ...        % T2s decay during adc
                              exp(-1i * ( phi_off*j + ...     % dephasing: offresonance 
                                          x * kx(:,j)' + ...  % dephasing: x-gradient
                                          y * ky(:,j)' ));    % dephasing: y-gradient
                DATA(j,:,:) = temp_signal.' * cmaps1D;
        end
        DATA = permute(DATA, [3,2,1]);

end

temp_t = toc;
fprintf(['   ' num2str(temp_t, '%.1f') 's ... complete! \n']);

end

