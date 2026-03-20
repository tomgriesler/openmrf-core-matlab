%% compare: sigpy SLR vs pulseq sinc
clear
pulseq_init();
FOV.dx = 1*1e-3;
FOV.dy = 1*1e-3;
FOV.dz = 5*1e-3;

INV.rf_type = 'HYPSEC_inversion';
INV.tExc    = 10 *1e-3;  % [s] pulse duration
INV.mu      = 4.9;         % [ ] determines amplitude of frequency sweep
INV.beta    = 600;       % [Hz] peak amplitude
INV = INV_init(INV, FOV, system);
INV_add();
seq.plot()

%% bloch simulation

df = -2000 : 50 : 2000;
n  = numel(df);
dt = system.rfRasterTime;

rf    = INV.rf.signal;
f1    = abs(rf);
phi   = phase(rf);
M_sim = zeros(n, 3);
parfor j=1:n    
    dw_ = 2*pi*df(j);
    M_  = [0; 0; 1];    
    for k=1:numel(rf)
        w1x_ = 2*pi * f1(k) * cos(phi(k));
        w1y_ = -2*pi * f1(k) * sin(phi(k));
        B_   = [  0     dw_    w1y_;
                 -dw_   0      w1x_;
                 -w1y_  -w1x_  0 ];
        M_ = expm(B_*dt) * M_;
    end    
    M_sim(j,:) = M_(:);
end
Mz  = M_sim(:,3);
Mxy = sqrt( M_sim(:,1).^2 + M_sim(:,2).^2 );
alpha = heaviside(sign(Mz)) .* asin(Mxy) + heaviside(-sign(Mz)) .* (pi-asin(Mxy));
% clear f1 phi M_sim Mz Mxy

%%
figure()
plot( df, alpha*180/pi, '-')
xlabel(['frequency [Hz]'])
ylabel(['flip angle [deg]'])
xline(0, 'b-')
yline(180,'k-')

