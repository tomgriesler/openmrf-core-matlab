%% calculate goia_wurst modulations
clear
pulseq_init();
system.maxGrad = 35 * 1e-3 * system.gamma;
clearvars -except system seq;

dur     = 3.5*1e-3;
f       = 0.9;
n_b1    = 16;
m_grad  = 4;
b1_max  = 650;
dz      = 5 *1e-3;
t_rise  = 2.5*1e-4;
t_spoil = 1*1e-3;


[rf, gz] = SIGPY_GOIA_WURST(dur, f, n_b1, m_grad, b1_max, dz, 0, t_rise, t_spoil, system);

seq.addBlock(rf, gz);
seq.plot();

%% interpolate rf and gz waveform for simulation
seq_rf = rf.signal.';
seq_gz = gz.waveform;
seq_gz = interp1(linspace(0,1,numel(seq_gz)), seq_gz, linspace(0,1,numel(seq_gz)*10))';
seq_rf = [seq_rf; zeros(numel(seq_gz)-numel(seq_rf),1)];
seq_rf = circshift(seq_rf, rf.delay/system.rfRasterTime);

seq_w1x = 2*pi * double(real(seq_rf));
seq_w1y = 2*pi * double(imag(seq_rf));
seq_gz  = 2*pi * double(seq_gz);

%% simulation of slice profile
Niso  = 10000;
z     = linspace(-1/2, 1/2, Niso)' * dz * 3;
R1    = 0;
R2    = 0;
Minit = [0; 0; 1; 1];

M = zeros(4, Niso);

db1 = 1;
df0 = 0;

parfor j=1:Niso
    M(:,j) = mex_BLOCH_rk4(Minit, seq_w1x*db1, seq_w1y*db1, z(j) * seq_gz + df0*2*pi, R1, R2, system.rfRasterTime);
end

nmean = 100;
Mx = M(1,:)';
My = M(2,:)';
Mz = M(3,:)';

figure()
hold on
plot(z*1e3, Mx, '-')
plot(z*1e3, My, '-')
plot(z*1e3, Mz, '-')
xline(-dz/2*1e3, 'k--', 'LineWidth', 2)
xline(dz/2*1e3,  'k--', 'LineWidth', 2)
