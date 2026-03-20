%%
clear
pulseq_init();
clearvars -except system seq;

FOV.Nxy      = 256;         % [ ] matrix size
FOV.Nz       = 1;           % [ ] numer of "stack-of-spirals", 1 -> 2D
FOV.fov_xy   = 256  *1e-3;  % [m] FOV geometry
FOV.dz       = 5   *1e-3;   % [m] slab or slice thickness
FOV.z_offset = 0    *1e-3;  % [m] slice offset
FOV.fov_z    = FOV.dz;
FOV_init();

ADIASL.mode  = 'on';
ADIASL.N_HS  = [1 2 4 6 8]';
ADIASL.tau   = 0.015;
ADIASL.f1max = 500;
ADIASL.beta  = 5;
ADIASL.dfmax = 200;
ADIASL       = ADIASL_init(ADIASL, FOV, system);

for loop_ADIASL = 1:numel(ADIASL.N_HS)
    ADIASL_add();
    seq.addBlock(mr.makeDelay(0.1));
end

seq.plot()

%% bloch simulation

% phase_cycle = [0]; % HS-1 (inversion)
% phase_cycle = [0 1]; % HS-2
% phase_cycle = [0 1 1 0]; % HS-4
% phase_cycle = [0 1 1 0   0 1]; % HS-6
phase_cycle = [0 1 1 0   1 0 0 1]; % HS-8

f1 = [];  
for temp_n_hs = 1:numel(phase_cycle)
    temp_rf = ADIASL.rf;
    if mod(temp_n_hs,2)==1        
        temp_rf.signal = ADIASL.f1 .* exp(1i*ADIASL.phi_down);
    else
        temp_rf.signal = ADIASL.f1 .* exp(1i*ADIASL.phi_up);
    end
    temp_rf.signal = temp_rf.signal * exp(1i*phase_cycle(temp_n_hs)*pi);        
    f1 = [f1; temp_rf.signal];    
end
clear temp_n_hs temp_rf;

df0 = -250 : 5 : 250;
db1 = 0.01 : 0.025 : 1.5;

n_df0 = numel(df0);
n_db1 = numel(db1);
Mz    = zeros(n_db1, n_df0);

parfor y=1:n_db1
for    x=1:n_df0
    M = mex_BLOCH_rk4( [0;0;1;1], db1(y)*2*pi*real(f1), db1(y)*2*pi*imag(f1), 2*pi*df0(x)*ones(size(f1)), 0, 0, system.rfRasterTime);
    Mz(y,x) = M(3);
end
end

figure()
imagesc(df0, (db1-1)*100, Mz, [-1 1]); axis square; colormap(jet); colorbar;
xlabel('df0 off-resonance [Hz]');
ylabel('db1 deviation [%]');
title(['Mz: HS-' num2str(numel(phase_cycle))])
set(gca, 'YDir', 'normal', 'FontName', 'arial', 'FontWeight', 'bold', 'FontSize', 12);