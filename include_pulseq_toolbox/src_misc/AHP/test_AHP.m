%% test adiabatic half passage (AHP)
clear
pulseq_init();

% old optimization for 3ms 600Hz
% AHP.tau   = 3*1e-3;
% AHP.wmax  = 600*2*pi;
% AHP.ratio = 0.960659060632200;
% AHP.beta1 = 0.436333936716347;
% AHP.beta2 = 0.787647050084116;

% new optimization
AHP = AHP_get_params(3*1e-3, 600*2*pi);

% get rf object
rf  = AHP_pulse( system, 'tipdown', AHP.tau, AHP.wmax, AHP.ratio, AHP.beta1, AHP.beta2, 0 );
% rf  = AHP_pulse( system, 'tipup',   AHP.tau, AHP.wmax, AHP.ratio, AHP.beta1, AHP.beta2, 0 );

seq.addBlock(rf);
seq.plot()

%% simulate

db1 = 0.75 : 0.05 : 1.25;
df0 = -50  : 10   : 50;
[db1, df0] = meshgrid(db1, df0);
db1 = db1(:);
df0 = df0(:);

w1x  = 2*pi*real(rf.signal);
w1y  = 2*pi*imag(rf.signal);
n_dt = numel(w1x);
n_sp = numel(db1);

M = zeros(3, n_sp, n_dt);

parfor j=1:n_sp
    M_ = [0; 0; 1];
    for k=1:n_dt
        w1x_ = w1x(k) * db1(j);
        w1y_ = w1y(k) * db1(j);
        dw0_ = 2*pi * df0(j);            
        B_   = [ 0,      dw0_,  -w1y_;
                -dw0_,   0,      w1x_;
                 w1y_,  -w1x_,   0];
        M_ = expm(B_*system.rfRasterTime) * M_;
        M(:,j,k) = M_;
    end
end

figure()
subplot(1,3,1)
plot(squeeze(M(1,:,:))')
subplot(1,3,2)
plot(squeeze(M(2,:,:))')
subplot(1,3,3)
plot(squeeze(M(3,:,:))')
