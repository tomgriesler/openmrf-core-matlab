%% test BIR4 pulse
clear
pulseq_init();
beta           = 10;
kappa          = atan(10);
dw0            = 30000;
f1             = 640;
alpha          = 90 *pi/180;
dphi           = 0.175;
pulse_duration = 10 *1e-3;
mod            = 'tipdown';
phase_offset   = 0;

[rf] = SIGPY_BIR4(beta, kappa, dw0, f1, alpha, dphi, pulse_duration, phase_offset, mod, system);

figure()
subplot(1,2,1)
plot(rf.t*1e3, abs(rf.signal))
xlabel('time t [ms]')
ylabel('abs(f1) [Hz]')

subplot(1,2,2)
plot(rf.t*1e3, phase(rf.signal))
xlabel('time t [ms]')
ylabel('angle(f1) [rad]')

f1  = rf.signal * exp(1i*phase_offset);
w1x = 2*pi * real(f1);
w1y = 2*pi * imag(f1);
dt  = system.rfRasterTime;

if strcmp(mod, 'tipdown')
    M_  = [0; 0; 1];
elseif strcmp(mod, 'tipup')
    M_  = [1; 0; 0];
else
    error('wrong mode: tipdown or tipup')
end
   
for j=1:numel(rf.signal)
    B_   = [  0        0      -w1y(j);
              0        0       w1x(j);
              w1y(j)  -w1x(j)  0 ];
    M_ = expm(B_*dt) * M_;
    M(:,j) = M_(:);
end

figure()
subplot(1,3,1)
plot3(M(1,:),M(2,:),M(3,:))
subplot(1,3,2)
plot(M(1,:),M(2,:))
xlabel('Mx')
ylabel('My')
subplot(1,3,3)
plot(rf.t*1e3, angle(f1)*180/pi)
xlabel('time t [ms]')
ylabel('angle [deg]')

disp(' ')
disp('final magnetization M:');
disp(['Mx: ' num2str(M_(1,end), '%.2f')]);
disp(['My: ' num2str(M_(2,end), '%.2f')]);
disp(['Mz: ' num2str(M_(3,end), '%.2f')]);
disp(' ')