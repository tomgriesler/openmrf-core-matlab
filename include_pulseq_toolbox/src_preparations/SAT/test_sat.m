%% compare: sigpy SLR vs pulseq sinc
clear
pulseq_init();
FOV.dx = 1*1e-3;
FOV.dy = 1*1e-3;
FOV.dz = 5*1e-3;

SAT.mode    = 'on';
SAT.rf_type = 'adiabatic_BIR4';
SAT = SAT_init(SAT, FOV, system);
SAT_add();
seq.plot()

%% bloch simulation

df0 = -100 : 10 : 100;
db1 = 0.75 : 0.025 : 1.25;
n_df0 = numel(df0);
n_db1 = numel(db1);
[df0, db1] = meshgrid(df0, db1);

dt  = system.rfRasterTime;
rf  = SAT.rf.signal; rf(abs(rf)==0) = [];
f1  = abs(rf);
phi = phase(rf);
Mz  = zeros(n_db1, n_df0);
Mxy = zeros(n_db1, n_df0);

parfor y = 1:n_db1
for    x = 1:n_df0
    dw_  = 2*pi*df0(x);
    db1_ = db1(y);
    M_   = [0; 0; 1];    
    for k=1:numel(rf)
        w1x_ = 2*pi * f1(k) * cos(phi(k)) * db1_;
        w1y_ = 2*pi * f1(k) * sin(phi(k)) * db1_;
        B_   = [  0     dw_    w1y_;
                 -dw_   0      w1x_;
                 -w1y_  -w1x_  0 ];
        M_ = expm(B_*dt) * M_;
    end
    Mxy(y,x) = sqrt(M_(1)^2+M_(2)^2);
    Mz(y,x)  = M_(3);
end
end

alpha = heaviside(sign(Mz)) .* asin(Mxy) + heaviside(-sign(Mz)) .* (pi-asin(Mxy));

%%
figure()
imagesc(alpha*180/pi, [0.75 1.25]*90); axis image; axis off; colormap jet; colorbar;

