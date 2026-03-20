% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% fixed parameters
acq.accfac    = 4;
acq.dninc     = 2;
acqP.flipex   = 90*pi/180;
acqP.MTfac    = 1;
acqP.nrep     = 1;
acqP.NSlices  = 1;
acqP.nTE      = 1;
acqP.sat_ppm  = -3.45;
acqP.sliceGAP = 1.25;
acqP.TEprep   = 10e-3;
acqP.TRfac    = [1 1 3 4 6 8];

seg.n_acc   = 1.4;
seg.n_full  = 1;
segP.fspS   = 0.5;
segP.GSfac  = 2;
segP.GXfac  = 1;
segP.GYfac  = 0;
segP.t1     = 1e-3;
segP.t2     = 0.5e-3;
segP.TEprep = 0.9e-3;
segP.tSp    = 1.5e-3;

spiral.kOffset = 50;
spiral.N       = 0.7;

dG     = 200e-6;
count  = 10;
kconc  = 1;
rf_fst = 8e-3;
SLfac  = 1;