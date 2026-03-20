%% test sigpy hyperbolic secans pulse
clear
pulseq_init();
beta           = 800;
mu             = 4.9;
pulse_duration = 10 *1e-3;
phase_offset   = 0;

[rf_sigpy] = SIGPY_HYPSEC(beta, mu, pulse_duration, phase_offset, system);

figure()
subplot(1,2,1)
plot(rf_sigpy.t*1e3, abs(rf_sigpy.signal))
xlabel('time t [ms]')
ylabel('abs(f1) [Hz]')

subplot(1,2,2)
plot(rf_sigpy.t*1e3, angle(rf_sigpy.signal)*180/pi)
xlabel('time t [ms]')
ylabel('angle(f1) [deg]')

%% bloch sim

df0 = 100 * linspace(-1, 1, 11);     % [Hz] offresonance
db1 = 0.2 * linspace(-1, 1, 11) + 1; % [ ]  B1+ error

Mz_final = zeros(numel(df0), numel(db1));

figure(345)
hold on

for l1 = 1:numel(df0)
for l2 = 1:numel(db1)

M      = [0; 0; 1];
Mx_sol = zeros(numel(rf_sigpy.signal),1);
My_sol = zeros(numel(rf_sigpy.signal),1);
Mz_sol = zeros(numel(rf_sigpy.signal),1);

for j = 1:numel(rf_sigpy.signal)
    dw0_ = 2*pi * df0(l1);
    w1x_ = 2*pi * real(rf_sigpy.signal(j)) *db1(l2);
    w1y_ = 2*pi * imag(rf_sigpy.signal(j)) *db1(l2);
    B_   = [  0      dw0_  -w1y_;
             -dw0_   0      w1x_;
              w1y_  -w1x_   0 ];
    M = expm(B_*system.rfRasterTime) * M;    
    Mx_sol(j) = M(1);
    My_sol(j) = M(2);
    Mz_sol(j) = M(3);

end

Mz_final(l1,l2) = M(3);

plot(rf_sigpy.t*1e3, Mz_sol, 'LineWidth', 1)
ylim([-1.1 1.1])

end
end

xlabel('time t [ms]')
ylabel('Mz [-1 ... 1]')
