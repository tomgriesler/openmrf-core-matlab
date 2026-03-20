clear
pulseq_init();
clearvars -except system

% input parameters
SPI.geo.design    = 'cones';
SPI.geo.N_loops   = 16;
SPI.geo.theta     = pi/8;
SPI.geo.FOV_coeff = [1 -0.5];
SPI.geo.kmax      = 256 / (2*0.256);
SPI.geo.Ns        = 1e5;
SPI.geo.ds        = 1e-3;
SPI.geo.flag_rv   = 0;

% limit read gradients
SPI.geo.lim_grad = 1.0;
SPI.geo.lim_slew = 1.0;

% calc geometry object
SPI.geo = SPI_minTimeCones(SPI.geo, system, 1);

%%
N = 48;
[sph, xyz, area] = SPI_fibonacci_sphere(N, 1);

figure('Color','w');
hold on
[Xs,Ys,Zs] = sphere(100);
surf(Xs,Ys,Zs,'EdgeColor','none','FaceAlpha',0.1);
for j=1:N
    az    = sph(j,1);
    colat = sph(j,2);   
    cph   = cos(az);
    sph_  = sin(az);
    cth   = cos(colat);
    sth   = sin(colat);    
    R     = [ cph*cth,  -sph_,  cph*sth;
              sph_*cth,  cph,   sph_*sth;
             -sth,       0,     cth ];
    krot = (R * SPI.geo.k')';
    plot3(krot(:,1), krot(:,2), krot(:,3))
end
axis equal off; view(35,20); camva('manual');
