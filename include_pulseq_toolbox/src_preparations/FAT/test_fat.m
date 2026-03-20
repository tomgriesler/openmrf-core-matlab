%% compare: sigpy SLR vs pulseq sinc
clear
pulseq_init();

FOV.dx = 5*1e-3;
FOV.dy = 5*1e-3;
FOV.dz = 5*1e-3;

system.B0 = 1.5;

FAT.mode    = 'on';
FAT.rf_type = 'slr';
FAT.tExc    = 8  *1e-3;
FAT.alpha   = 100 *pi/180;
FAT = FAT_init(FAT, FOV, system);
FAT_add();
seq.plot()

%% bloch simulation

df = -2*abs(FAT.rf.freqOffset) : 10 : 2*abs(FAT.rf.freqOffset);
n  = numel(df);
dt = system.rfRasterTime;

rf    = FAT.rf.signal;
f1    = abs(rf);
phi   = phase(rf);
M_sim = zeros(n, 3);
parfor j=1:n    
    dw_ = 2*pi*df(j);
    M_  = [0; 0; 1];    
    for k=1:numel(rf)
        w1x_ = 2*pi * real(rf(k));
        w1y_ = 2*pi * imag(rf(k));
        B_   = [  0     dw_    -w1y_;
                 -dw_   0      w1x_;
                 w1y_   -w1x_  0 ];
        M_ = expm(B_*dt) * M_;
    end    
    M_sim(j,:) = M_(:);
end
Mz  = M_sim(:,3);
Mxy = sqrt( M_sim(:,1).^2 + M_sim(:,2).^2 );
alpha = heaviside(sign(Mz)) .* asin(Mxy) + heaviside(-sign(Mz)) .* (pi-asin(Mxy));
clear f1 phi M_sim Mz Mxy

%%
figure()
plot( df + abs(FAT.rf.freqOffset), alpha*180/pi, '-')
xlabel(['frequency [Hz]'])
ylabel(['flip angle [deg]'])
xline(0, 'b-')
xline(abs(FAT.rf.freqOffset)/2, 'r--')
xline(abs(3*FAT.rf.freqOffset)/2, 'r--')
yline(90,'k-')
ylim([0, 110])
