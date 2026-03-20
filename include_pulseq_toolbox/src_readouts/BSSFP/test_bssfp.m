clear

pulseq_init();

%% FOV geometry
FOV.Nx       = 256; % matrix sixe: x
FOV.Ny       = 192; % matrix sixe: y
FOV.Nz       = 96;  % matrix sixe: z
FOV.dz       = 128 *1e-3; % [m] slice/slab thickness
FOV.fov_x    = 256 *1e-3; % [m] FOV: x direction
FOV.fov_y    = 256 *1e-3; % [m] FOV: y direction
FOV.fov_z    = 128 *1e-3; % [m] FOV: z direction
FOV.z_offset = 0 *1e-3;   % [m] shift in z direction (can be done at the scanner, leave 0)
FOV_init();

%% BSSFP sequence parameters
BSSFP.TR        = 20  *1e-3;  % [s] repetition time, use 0 for auto-minimization
BSSFP.t_acq     = 3.2 *1e-3;  % [s] acquisition time
BSSFP.t_pre_rew = 2.0 *1e-3;  % [s] duration of prephaser and rewinder gradients
BSSFP.os_mode   = 0;          % [0/1] read oversampling: 0 off, 1 on

% slice excitation
BSSFP.exc_mode       = 'sinc';      % 'sinc' or 'sigpy_SLR'
BSSFP.exc_shape      = 'ex';        % 'st' or 'ex' only for sigpy_SLR case 
BSSFP.exc_time       = 2.0 *1e-3;   % [s] excitation time
BSSFP.exc_tbw        = 6;           % [ ] time bandwidth product 
BSSFP.exc_flipangle  = 20 *pi/180;  % [rad] flip angles

BSSFP = BSSFP_init(BSSFP, FOV, system);

%% write example sequence
seq = mr.Sequence(system);
for loop_z = 1:FOV.Nz
    for loop_y = 1:FOV.Ny
        BSSFP_add();
    end
end
