function [seq, SPITSE] = SPITSE_add(seq, system, FOV, SPITSE, loop_dummy)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% Single shot spiral TSE with annulated segmentation
% https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.29224

    % The |Sequence| class provides functionality to create magnetic
    % resonance sequences (MRI or NMR) from basic building blocks.
    %
    % This provides an implementation of the open file format for MR sequences
    % described here: http://pulseq.github.io/specification.pdf
    %
    % This example performs the following steps:
    %
    % # Create slice selective RF pulse for imaging
    % # Create readout gradient and phase encode strategy.
    % # Loop through phase encoding and generate sequence blocks.
    % # Write the sequence to an open file format suitable for execution on a
    % scanner.
    %
    % The sequence runs in different modes defined by :
    %   scanmode=
    %       'init'  : initialization, creates empty k-space corrections dKA and dKE
    %       'trim'  : used for experimental trimming,
    %                 the first and every 2:decho:acqP.necho echo periods are
    %                 readout as a spin echo
    %       'run'   : run sequence
    %
    %   segmode=
    %       'fix'   : fixed number of points per segment
    %       'tan'   : tangential transition to and from each segment
    %       'single': single echo spiral
    %
    %   spmode=
    %       'cont'  : continuous segmentation
    %       'alt'   : alternate between spiral-out and spiral-in in odd and
    %                 even segments.
    %
    %   initmode=
    %       'no'    : no gradient in first half of first segment
    %       'rev'   : first half is reverse of second half
    %       'outin' : first half 180ï¿½ rotated and out-in
    %       'sym'   : first half mirrored at y-axis and out-in;
    %
    %   accmode=    acceleration mode
    %       'dd'    : dual density. First seg.n_full*acq.nadc/2 points are fully sampled,
    %                 others undersampled by a factor seg.n_acc
    %       'vd'    : variable density with Fcoeff2 = - Fcoeff/seg.n_acc
    %
    %   fatsat
    %       'on'
    %       'no'    :  with fatsat gradient, but no rf
    %       'off'
    %   T2prep
    %       'on'    : First echo is different from CPMG-train. Usefúll for fMRI
    %                and quantitative T2-measurements
    %
    %   seq_var_mod     switch for specific parameters. Some of the switches require to set kseq>1
    %                   Parameters for the different modes are defined in mySpiralTSE_par.m
    %                   
    %       seqvar_mod='none';      normal execution
    %       seqvar_mod='T1var';     generates sequence with multiple
    %                               acquisition of individual slices for T1-calculation
    %                               acquisition scheme is determined by
    %                               acqP.TRfac, which is an array representing
    %                               the slice interleaving scheme.                            
    %       seqvar_mod='MTvar';     used for assessment of MT-effect. acqP.MTfac determines the position of central slice       
    %       seqvar_mod='TErep';     number of images acquired on a single
    %       echotrain.              Used to measure long T2.
    %       seqvar_mod='TEvar';     Generates sequence with TE-variation by concatenation
    %       
    %   plotflag:   sequence of 6 binary digits to define, which plots are shown
    %       1   :   plots from vds-package
    %       2   :   trajectory as cloud of dots to check completeness
    %       3   :   plot of individual segments
    %       4   :   x-y-trajectory during readout(red) and total trajectory (blue)
    %       5   :   x-z-trajectory during readout(red) and total trajectory (blue)
    %       6   :   sequence

% only for calc trajectories
if isempty(seq)
    seq = mr.Sequence(system);
end

% dummy mode
if nargin<5
    loop_dummy = 1; % 0 -> dummy   1 -> adc
end

B0    = system.B0;
gamma = system.gamma;

%% load sequence parameters from SPITSE struct

% default values
if ~isfield(SPITSE, 'scanmode')
    SPITSE.scanmode = 'run';
end
if ~isfield(SPITSE, 'segmode')
    SPITSE.segmode = 'fix';
end
if ~isfield(SPITSE, 'spmode')
    SPITSE.spmode = 'cont';
end
if ~isfield(SPITSE, 'initmode')
    SPITSE.initmode = 'no';
end
if ~isfield(SPITSE, 'accmode')
    SPITSE.accmode = 'vd';
end
if ~isfield(SPITSE, 'fatsat')
    SPITSE.fatsat = 'off';
end
if ~isfield(SPITSE, 'seqvar_mod')
    SPITSE.seqvar_mod = 'none';
end
if ~isfield(SPITSE, 'T2prep')
    SPITSE.T2prep = 'on';
end
if ~isfield(SPITSE, 'EncMode')
    SPITSE.EncMode = 'pd_linear';
end
if ~isfield(SPITSE, 'plotflag')
    SPITSE.plotflag = '000';
end
if ~isfield(SPITSE, 'dispflag')
    SPITSE.dispflag = 0;
end
if ~isfield(SPITSE, 'slewfac')
    SPITSE.slewfac = 0.99;
end
if ~isfield(SPITSE, 'gradfac')
    SPITSE.gradfac = 0.99;
end
if ~isfield(SPITSE, 'tbw')
    SPITSE.tbw = 4;
end
if ~isfield(SPITSE, 'TR')
    SPITSE.TR = system.gradRasterTime;
end

%
scanmode   = SPITSE.scanmode;
segmode    = SPITSE.segmode;
spmode     = SPITSE.spmode;
initmode   = SPITSE.initmode;
accmode    = SPITSE.accmode;
fatsat     = SPITSE.fatsat;
seqvar_mod = SPITSE.seqvar_mod;
T2prep     = SPITSE.T2prep;
plotflag   = SPITSE.plotflag;
dispflag   = SPITSE.dispflag;

% load params from FOV struct
acqP.fov            = FOV.fov_xy;
spiral.Nx           = ceil(FOV.Nxy / 2);
acqP.sliceThickness = FOV.dz;
acqP.sliceOffset    = FOV.z_offset;

% load params from SPITSE struct
acqP.TR       = SPITSE.TR;
acqP.TE       = SPITSE.TE;
acqP.necho    = SPITSE.NEcho;
rfex_phase    = SPITSE.rfex_phase;
rfref_phase   = rfex_phase-pi/2;
acqP.PEtype   = SPITSE.EncMode;
segP.tEx      = SPITSE.tEX;
segP.tRef     = SPITSE.tRef;
acqP.flipref  = SPITSE.flipref; 
acqP.flipflag = SPITSE.flipflag;
acq.slewfac   = SPITSE.slewfac;
acq.gradfac   = SPITSE.gradfac;

%% fixed parameters
SPITSE_params();

%% ---------- calculate spiral turbo spin echo blocks ----------
seqname = 'TSE';
warning('OFF', 'mr:restoreShape')

for loop_kseq = 1:1    
    count = count+1;    
    acqP.sat_freq = acqP.sat_ppm*1e-6*B0*system.gamma;    
    temp = sprintf('%s %s','fatsat:',fatsat);
    if dispflag==1
        disp(temp);
    end
    temp = sprintf('%s %s','T2prep:',T2prep);
    if dispflag==1
        disp(temp);
    end
    if(strcmp(scanmode,'trim'))
        seqname = strcat(seqname,'_',segmode,'_',spmode,'_trim');
    end
    acqP.sliceThickness = 4e-3*SLfac;
    acqP.sliceGAP = acqP.sliceGAP/SLfac;    
    acqP.sltemp = [1:2:acqP.NSlices 2:2:acqP.NSlices];
    if(strcmp(seqvar_mod,'T1var'))
        [temp]  =  T1var_index(acqP.NSlices,acqP.TRfac);
        acqP.sltemp = temp;
        acqP.Oslices = acqP.sltemp-ceil(acqP.NSlices/2);
        acqP.NSlices = length(temp);
    else
        if(strcmp(seqvar_mod,'MTvar'))
            [temp]  =  MTvar_index(acqP.sltemp,acqP.MTfac);
            acqP.sltemp = temp;
            acqP.Oslices = acqP.sltemp-ceil(acqP.NSlices/2);
            acqP.NSlices = length(temp);
        else
            acqP.Oslices = acqP.sltemp-ceil(acqP.NSlices/2);
        end
    end
    if(strcmp(T2prep,'off'))
        acqP.TEeff = acqP.nTE*acqP.TE;
    else
        acqP.TEeff = acqP.TEprep;
        acqP.nTE = 1;
    end
    segP.tExwd = segP.tEx+system.rfRingdownTime+system.rfDeadTime;
    segP.tRefwd = segP.tRef+system.rfRingdownTime+system.rfDeadTime;
    segP.tSpex = 0.5*(acqP.TE-segP.tExwd-segP.tRefwd);               
    n_acc_s = num2str(seg.n_acc);
    if(seg.n_acc == floor(seg.n_acc))
        n_acc_s = num2str(seg.n_acc);
    else
        ind = find(n_acc_s == '.'); n_acc_s(ind) = '_';
    end
    seqtype = strcat(sprintf('%03s',num2str(count)),'TSE');
    if(kconc == 2)
        seqtype = strcat(sprintf('%03s',num2str(count)),'TSC');
    end    
    if(loop_kseq>0)
        quitflag = 0;
        switch seqvar_mod
            case 'none'
                seqname = strcat('TSE_',T2prep,'_',num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices));
            case 'Diff'
                Diffvals = [4 10 20 30 40 4 10 20 30 40];
                if(loop_kseq)>10
                    quitflag = 1;
                    return;
                end
                segP.GXfac = Diffvals(loop_kseq);
                segP.GYfac = 0;
                if(loop_kseq>5) 
                    segP.GYfac = Diffvals(loop_kseq);
                    segP.GXfac = 0;
                end
                seqname = strcat('TSE_',T2prep,'_',num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices),'_Diff');
                allname = strcat('TSEall_',num2str(round(1000*acqP.TE)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices),'_Diff'); 
            case 'TErep'
                seqname = strcat('TSE_',T2prep,'_',num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices),'_',num2str(acqP.nrep));
            case 'T1var'
                seqname = strcat('TSE_T1var_',T2prep,'_',num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices));
            case 'MTvar'
                seqname = strcat('TSE_MTvar_',T2prep,'_',num2str(acqP.MTfac),'_',num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices));
            case 'TEvar'
                TEvals = [1 1 2 3 4 7 10 13 16];
                if(loop_kseq)>14
                    quitflag = 1;
                    return;
                end
                acqP.TEeff = TEvals(loop_kseq)*acqP.TE;
                acqP.TEprep = acqP.TEeff;            
                seqname = strcat('TSE_',T2prep,'_',num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices));
                allname = strcat('TSEvar_',num2str(round(1000*acqP.TE)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices)); 
            case 'TEall'
                TEvals = [40 40 60 80 100 140 [10:10:220]];
                if(loop_kseq)>28
                    quitflag = 1;
                    return;
                end
                acqP.TEeff = TEvals(loop_kseq)*10^-3;
                acqP.TEprep = acqP.TEeff;
                if(loop_kseq<6)
                    T2prep = 'on';
                    acqP.nTE = 1;
                else
                    T2prep = 'off';
                    acqP.nTE = round(acqP.TEeff/acqP.TE);
                end
                seqname = strcat('TSE_',T2prep,'_',num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices));
                allname = strcat('TSEall_',num2str(round(1000*acqP.TE)),'_',num2str(100*acqP.sliceGAP),'_',num2str(acqP.NSlices)); 
            case 'accvar'
                if(loop_kseq)>8, quitflag = 1;
                    return;
                end
                N_var = [0.7 1];
                acc_var = [1.2 1.4 2 4];
                if(2*floor(loop_kseq/2) == loop_kseq)
                    spiral.N = N_var(2); 
                else
                    spiral.N = N_var(1); 
                end
                seg.n_acc = acc_var(floor(loop_kseq/2));
                seqname = strcat(seqtype,num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(acqP.flipref),'_',num2str(10*seg.n_acc),'_',num2str(10*spiral.N),'_',num2str(acqP.NSlices));
            case 'fsvar'
                fatsat = 'on';
                if(loop_kseq)>7
                    quitflag = 1;
                    return;
                end
                fac = [1 1 1.1 1.2 1.3 1.4 1.5];
                acqP.sat_freq = fac(loop_kseq)*acqP.sat_ppm*1e-6*B0*system.gamma;
                seqname = strcat(seqtype,num2str(round(1000*acqP.TE)),'_',num2str(floor(1000*acqP.TEeff)),'_',num2str(acqP.flipref),'_',num2str(floor(abs(acqP.sat_freq))),'_',num2str(acqP.NSlices));
        end
        if(quitflag == 1)
            return;
        end
    end
    if(exist(strcat(seqname,'.seq'),'file')>0)
        m = 1; testname = seqname;
        while(exist(strcat(testname,'.seq'),'file')>0)
            testname = strcat(seqname,'_',num2str(m));
            m = m+1;
        end
        seqname = testname;
    end
    acq.readoutTime  =  acqP.TE-2*segP.tSp-segP.tRefwd-2*20e-6;
    [rfex, gz]  =  mr.makeSincPulse( acqP.flipex, system, 'Duration', segP.tEx, ...
                                     'sliceThickness', acqP.sliceThickness, ...
                                     'apodization', 0.5, ...
                                     'timeBwProduct', SPITSE.tbw, ...
									 'use', 'excitation', ...
                                     'maxSlew', system.maxSlew*4);
    GSex  =  mr.makeTrapezoid('z',system,'amplitude',gz.amplitude,'FlatTime',segP.tExwd,'riseTime',dG);
    [rfref, gz]  =  mr.makeSincPulse( pi, system, 'Duration', segP.tRef, ...
                                      'sliceThickness', acqP.sliceThickness, ...
                                      'apodization', 0.5, ...
                                      'timeBwProduct', SPITSE.tbw, ...
                                      'use', 'refocusing', ...
                                      'maxSlew', system.maxSlew*4);
    GSref  =  mr.makeTrapezoid('z',system,'amplitude',GSex.amplitude,'FlatTime',segP.tRefwd,'riseTime',dG);
    refenvelope = rfref.signal;
    AGSex = GSex.area/2;
    GSspr  =  mr.makeTrapezoid('z',system,'area',AGSex*(1+segP.fspS),'duration',segP.tSp,'riseTime',dG);
    GSspex  =  mr.makeTrapezoid('z',system,'area',AGSex*segP.fspS,'duration',segP.tSpex,'riseTime',dG);
    if (B0<2)
        rf_fst = 1e-5*floor(1e5*10e-3*1.5/B0);
    end
    rf_fs  =  mr.makeGaussPulse(110*pi/180,'system',system,'Duration',rf_fst,...
        'bandwidth',abs(acqP.sat_freq),'freqOffset',acqP.sat_freq,'use','saturation');
    gz_fs  =  mr.makeTrapezoid('z',system,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm
    spiral.deltak = 1/acqP.fov;
    spiral.kWidth  =  2*spiral.Nx*spiral.deltak;
    acq.ntot = round((acq.readoutTime+2*20e-6)/system.gradRasterTime);
    acq.nadc = acq.accfac*round((acq.readoutTime)/system.gradRasterTime);
    adc  =  mr.makeAdc(acq.nadc,'Duration',acq.readoutTime, 'Delay', 20e-6);%,'Delay',GRacq.riseTime);   
    smax  =  acq.slewfac*100*system.maxSlew/gamma;	 % 150 T/m/s
    gmax  =  acq.gradfac*100*system.maxGrad/gamma;	 % G/cm    
    T  =  system.gradRasterTime;	 % Seconds    
    Fcoeff  =  acqP.fov*100 ;
    if(strcmp(accmode,'vd'))
        Fcoeff = [acqP.fov*100 -acqP.fov*100/(seg.n_acc)];
    end    
    res = 100*acqP.fov/spiral.Nx;
    rmax  =  2*spiral.kWidth/100;		% cm^(-1), corresponds to 1mm resolution.
    rmax = 1/res;    
    if dispflag==1
        disp('Calculating Gradient');
    end
    [k,~,~,time,r,theta]  =  vds(smax,gmax,T/acq.accfac,spiral.N,Fcoeff,rmax);
    if ((spiral.N>1)&&strcmp(accmode,'dd'))
        [k1,~,~,time,r,theta]  =  vds(smax,gmax,T/acq.accfac,1,Fcoeff,rmax);
        k1 = k1(1:floor(seg.n_full*acq.nadc/2));
        k2 = k;
        rads = abs(k1(end));
        ind = find(abs(k2) >= rads);
        ang1 = angle(k1(end));
        ang2 = angle(k2(ind(1)));
        k2r = k2(ind(1)+1:end)*exp(i*(ang1-ang2));
        k = [k1 k2r];
    end    
    k = k(1:acq.accfac:end);    
    g  =  10^4/gamma*([k 0]-[0 k])/T;
    g  =  g(1:length(k));
    if dispflag==1
        disp('Plotting Gradient');
    end
    g  =  [real(g(:)) imag(g(:))];
    figure
    vds_plotgradinfo(g,T);
    if(plotflag(1) == '0')
        close;
    end
    np = length(k);
    kt = zeros([np 2]);
    kt(:,1) = real(k*100);
    kt(:,2) = imag(k*100);
    figure
    plot(kt(:,1),kt(:,2),'b.');
    if(plotflag(2) == '0')
        close;
    end
    Gsp = g*gamma/100;      
    acq.ninc = acq.ntot-2*acq.dninc;
    spiral.kStart = floor(spiral.kWidth/2)+spiral.kOffset;
    if(strcmp(segmode,'single'))
        acqP.necho = 1; seg.i1(1) = floor(acq.ntot/2);
        seg.ikseg = zeros([acqP.necho 2]);
        seg.ikseg(1,1) = 1; seg.ikseg(1,2) = floor(acq.ninc/2);
        seg.nkseg = acq.ninc*ones([acqP.necho 1]); seg.nkseg(1) = floor(acq.ninc/2);
    end
    if(strcmp(segmode,'tan'))
        [seg.kseg, seg.nkseg, seg.ikseg, acqP.necho, seg.i1]  =  spiral_seg(kt,acq.ninc,spiral.kStart,acq.ntot,spmode);
        axlim = max(spiral.kStart,spiral.kWidth/2);
        axis([-axlim axlim -axlim axlim])
        if(plotflag(3) == '0'), close; end
    end    
    if(strcmp(segmode,'fix'))        
        acqP.necho = floor(np/acq.ninc);
        seg.nkseg = acq.ninc*ones([acqP.necho 1]); seg.nkseg(1) = floor(acq.ninc/2);
        seg.i1 = 0*seg.nkseg+acq.dninc; seg.i1(1) = floor(acq.ntot/2);
        seg.ikseg = zeros([acqP.necho 2]);
        seg.ikseg(1,1) = 1; seg.ikseg(1,2) = floor(acq.ninc/2);        
        for k = 2:acqP.necho
            seg.ikseg(k,1) = seg.ikseg(k-1,2);
            seg.ikseg(k,2) = seg.ikseg(k,1)+seg.nkseg(k)-1;
        end
        if(strcmp(spmode,'alt'))
            temp = fliplr(seg.ikseg);
            seg.ikseg(2:2:end,:) = temp(2:2:end,:);
        end        
        if(plotflag(3) == '1')
            figure
            axlim = max(spiral.kStart,spiral.kWidth/2);
            axis([-axlim axlim -axlim axlim])            
            axis square
            hold on
            col = ['b-';'r-'];
            if(strcmp(spmode,'cont'))
                for k = 1:2:acqP.necho
                    plot(kt(seg.ikseg(k,1):seg.ikseg(k,2),1),(kt(seg.ikseg(k,1):seg.ikseg(k,2),2)),col(1,:));
                    plot(kt(seg.ikseg(k,2),1),kt(seg.ikseg(k,2),2),'bo');
                end
                col = flipud(col);
                for k = 2:2:acqP.necho
                    plot(kt(seg.ikseg(k,1):seg.ikseg(k,2),1),(kt(seg.ikseg(k,1):seg.ikseg(k,2),2)),col(1,:));
                    plot(kt(seg.ikseg(k,1),1),kt(seg.ikseg(k,1),2),'ro');
                    
                end
            end
            if(strcmp(spmode,'alt'))
                for k = 1:2:acqP.necho
                    plot(kt(seg.ikseg(k,1):seg.ikseg(k,2),1),(kt(seg.ikseg(k,1):seg.ikseg(k,2),2)),col(1,:));
                    plot(kt(seg.ikseg(k,2),1),kt(seg.ikseg(k,2),2),'bo');
                end
                col = flipud(col);
                for k = 2:2:acqP.necho
                    plot(kt(seg.ikseg(k,1):-1:seg.ikseg(k,2),1),(kt(seg.ikseg(k,1):-1:seg.ikseg(k,2),2)),col(1,:));
                    plot(kt(seg.ikseg(k,1),1),kt(seg.ikseg(k,1),2),'ro');                    
                end
            end
        end
    end           
    dKA = zeros([acqP.necho 2]);
    dKE = zeros([acqP.necho 2]);    
    
    % split gradients and recombine into blocks
    % lets start with slice selection....
    GS1times = [0 GSex.riseTime];
    GS1amp = [0 GSex.amplitude];
    GS1  =  mr.makeExtendedTrapezoid('z','times',GS1times,'amplitudes',GS1amp);
    
    GS2times = [0 GSex.flatTime];
    GS2amp = [GSex.amplitude GSex.amplitude];
    GS2  =  mr.makeExtendedTrapezoid('z','times',GS2times,'amplitudes',GS2amp);
    
    GS3times = [0 GSspex.riseTime GSspex.riseTime+GSspex.flatTime GSspex.riseTime+GSspex.flatTime+GSspex.fallTime];
    GS3amp = [GSex.amplitude GSspex.amplitude GSspex.amplitude GSref.amplitude];
    GS3  =  mr.makeExtendedTrapezoid('z','times',GS3times,'amplitudes',GS3amp);
    GS3p1 = mr.makeExtendedTrapezoid('z','times',[0 GSex.fallTime],'amplitudes',[GSex.amplitude 0]);
    
    %GS3p3 = mr.makeExtendedTrapezoid('z','times',[0 GSref.fallTime],'amplitudes',[GSref.amplitude 0]);
    GS4times = [0 GSref.flatTime];
    GS4amp = [GSref.amplitude GSref.amplitude];
    GS4  =  mr.makeExtendedTrapezoid('z','times',GS4times,'amplitudes',GS4amp);
    GS5times = [0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
    GS5amp = [GSref.amplitude GSspr.amplitude GSspr.amplitude 0];
    GS5  =  mr.makeExtendedTrapezoid('z','times',GS5times,'amplitudes',GS5amp);
    GS5p1times = [0 dG];
    GS5p1amp = [GSref.amplitude 0];
    GS5p1  =  mr.makeExtendedTrapezoid('z','times',GS5p1times,'amplitudes',GS5p1amp);

    GS7times = [0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
    GS7amp = [0 GSspr.amplitude GSspr.amplitude GSref.amplitude];
    GS7  =  mr.makeExtendedTrapezoid('z','times',GS7times,'amplitudes',GS7amp);
    GS7p1times = [0 GSspr.riseTime GSspr.riseTime+segP.TEprep GSspr.riseTime+segP.TEprep+GSspr.fallTime];
    GS7p1amp = [0 segP.GSfac/2*GSspr.amplitude segP.GSfac/2*GSspr.amplitude GSref.amplitude];
    GS7p1 =  mr.makeExtendedTrapezoid('z','times',GS7p1times,'amplitudes',GS7p1amp);
    
    GSarea = calcArea(GS4)/2+calcArea(GS7);
    GSprep_area1 = calcArea(GS2)/2+calcArea(GS3p1)+calcArea(GS7p1)+calcArea(GS4)/2;
    GSprep_area2 = calcArea(GS4)/2+calcArea(GS5p1);
    GS3p2 = mr.makeTrapezoid('z',system,'area',segP.GSfac*GSarea-GSprep_area1,'duration',segP.t1,'riseTime',dG);
    GS5p2 = mr.makeTrapezoid('z',system,'area',segP.GSfac*GSarea-GSprep_area2,'duration',(acqP.TE-2*segP.tSp-segP.tRefwd)/2,'riseTime',dG);
    
    % and now the readout gradient....
    GRpre_s  =  mr.makeTrapezoid('x',system,'area',spiral.kStart,'duration',segP.t1,'riseTime',dG);
    GRpre  =  mr.makeTrapezoid('x',system,'area',spiral.kStart,'duration',segP.tSpex,'riseTime',dG);
    GRref  =  mr.makeTrapezoid('x',system,'area',spiral.kStart,'duration',acq.readoutTime,'riseTime',dG);
    GRtrim  =  mr.makeTrapezoid('x',system,'area',2*spiral.kStart,'duration',acqP.TE-GSref.flatTime-2*segP.tSp,'riseTime',dG);
    GRspoi_x = mr.makeTrapezoid('x',system,'area',(segP.GXfac-1)*spiral.kStart,'duration',GSspr.riseTime+segP.TEprep,'riseTime',dG);
    GRspoi_xs = mr.makeTrapezoid('x',system,'area',(segP.GXfac-1)*spiral.kStart,'duration',GSspr.riseTime+segP.TEprep,'riseTime',dG);
    GRspoi_y = mr.makeTrapezoid('y',system,'area',segP.GYfac*spiral.kStart,'duration',GSspr.riseTime+segP.TEprep,'riseTime',dG);
    
    % and filltimes
    segS.tEx = GS1.shape_dur+GS2.shape_dur+GS3.shape_dur;
    segS.tRef = GS4.shape_dur+GS5.shape_dur+GS7.shape_dur+acq.readoutTime;
    tend = GS4.shape_dur+GS5.shape_dur;
    tETrain = segS.tEx+acqP.necho*segS.tRef+tend;
    TRfill = acqP.TR;
    if TRfill<0, TRfill = 1e-3;
        disp(strcat('Warning!!! acqP.TR too short, adapted to include all slices to : ',num2str(1000*acqP.NSlices*(tETrain+TRfill)),' ms'));
    else
        if dispflag==1
            disp(strcat('TRfill : ',num2str(1000*TRfill),' ms'));
        end
    end
    delayTR  =  mr.makeDelay(TRfill);
    if(strcmp(T2prep,'on'))
        TEfill1 = acqP.TEprep/2-GS2.shape_dur/2-GS3p1.shape_dur-segP.t1-GS7p1.shape_dur-GS4.shape_dur/2;
        TEfill2 = segP.t2+acqP.TEprep/2+acqP.TE/2-(segP.tRef+GS5p1.shape_dur+segP.tSp+acqP.TE-2*segP.tSp-segP.tRefwd+GS7.shape_dur);
        delayTE1  =  mr.makeDelay(TEfill1);
        delayTE2  =  mr.makeDelay(TEfill2);
    end
    % and flip angles
    acqP.rflip = acqP.flipref+zeros([1 acqP.necho]);
    if(acqP.flipflag == 1)
        acqP.rflip(1) = 90+acqP.flipref/2;
    end
    if(acqP.flipflag == 2)
        [rf,~]  =  fliptraps(acqP.flipref,(acqP.nrep+1)*acqP.necho,6,'opt',0,2,0,acqP.flipref,acqP.flipref,[6 5 5 acqP.necho]);
        rf = [180 rf];
        rf = rf(1:acqP.nrep*acqP.necho); pow = sum(rf.^2)/sum((0*rf+180).^2);
        if dispflag==1
            disp(strcat('rel. power :',num2str(pow)));
        end
        acqP.rflip = rf(1:acqP.nrep*acqP.necho);
    end
    if(strcmp(acqP.PEtype,'linear'))
        acqP.PEind = acqP.necho-(mod(acqP.nTE-1+(acqP.necho-[1:acqP.necho]),acqP.necho));
    end
    if(strcmp(acqP.PEtype,'centric'))
        acqP.PEind = zeros([1 acqP.necho]);
        acqP.PEind(1:acqP.nTE-1) = [2*acqP.nTE-2:-2:2];
        acqP.PEind(acqP.nTE:2*acqP.nTE-1) = [1:2:2*acqP.nTE-1];
        acqP.PEind(2*acqP.nTE:end) = [2*acqP.nTE:acqP.necho];
    end
    if(strcmp(acqP.PEtype,'pd_linear'))
        acqP.PEind = 1 : 1 : acqP.necho;
    end
    
    % Define sequence blocks
    if(~exist('scanmode'))
        scanmode = 'run';
    end
    ntr = 0;
    switch scanmode
        case 'init'
            dKA = zeros([acqP.necho 2]);
            dKE = zeros([acqP.necho 2]); % dKA and dKE are trim values to correct for trajectory imperfections
        case 'trim'
            decho = 4;
            acqP.sliceGAP = 0;
            ntrim = [1 2:decho:acqP.necho];
            acqP.NSlices = length(ntrim);
        case 'allTE'
            acqP.sliceGAP = 0;
            acqP.NSlices = acqP.necho;
    end
    nk = acqP.necho;    
    for loop_s = 1:acqP.NSlices
        dz                =                   acqP.sliceGAP * acqP.sliceThickness * (acqP.Oslices(loop_s));
        rfex.freqOffset   = GSex.amplitude  * acqP.sliceGAP * acqP.sliceThickness * (acqP.Oslices(loop_s)) + GSex.amplitude  * acqP.sliceOffset;
        rfref.freqOffset  = GSref.amplitude * acqP.sliceGAP * acqP.sliceThickness * (acqP.Oslices(loop_s)) + GSex.amplitude  * acqP.sliceOffset;
        rfex.phaseOffset  = rfex_phase-2*pi*rfex.freqOffset*mr.calcRfCenter(rfex); % align the phase for off-center slices
        rfref.phaseOffset = rfref_phase-2*pi*rfref.freqOffset*mr.calcRfCenter(rfref); % dito
        rfref.signal = refenvelope;            
        if(strcmp(fatsat,'on'))
            seq.addBlock(rf_fs,gz_fs);
        end
        if(strcmp(fatsat,'no'))
            seq.addBlock(gz_fs);
        end            
        if(strcmp(T2prep,'on'))
            seq.addBlock(GS1);
            seq.addBlock(GS2,rfex);
            seq.addBlock(GS3p1);
            seq.addBlock(GS3p2,GRpre_s);
            seq.addBlock(delayTE1);
            seq.addBlock(GS7p1,GRspoi_x,GRspoi_y);
            seq.addBlock(GS4,rfref);
            seq.addBlock(GS5p1);
            seq.addBlock(delayTE2);
        else
            seq.addBlock(GS1);
            seq.addBlock(GS2,rfex);
            seq.addBlock(GS3,GRpre);
        end
        if(ntr>0)
            nk = loop_s;
        end
        for loop_m = 1:acqP.nrep*acqP.necho
            if(strcmp(scanmode,'allTE'))
                acqP.nTE = loop_s;
            end
            indk = mod(loop_m,acqP.necho);
            if(indk == 0)
                indk = acqP.necho;
            end
            k = acqP.PEind(indk);
            if((k == 1)&&strcmp(initmode,'rev'))
                ind = find(kt(1:seg.nkseg(1),2)>0);
                nkseg1 = max(ind);
                segS.tSp = GSspr.riseTime+GSspr.flatTime+GSspr.fallTime+(seg.i1(k)-nkseg1)*system.gradRasterTime;
                kBegin1 = [-spiral.kStart+dKA(k,1) 0];
                kEnd1 = [kt(nkseg1,1) kt(nkseg1,2)+dKA(k,2)];
                GStart = [0 0];
                GDest = [-Gsp(seg.ikseg(k,2),1) -Gsp(seg.ikseg(k,2),2)];
                if(strcmp(segmode,'tan'))
                    [GBegin, tBegin, Ginit, tinit]  =  spiral_k2k_m(kBegin1, kEnd1, GStart, GDest, segS.tSp, system.maxSlew, system.gradRasterTime);
                else
                    [GBegin, tBegin, tt]  =  spiral_k2k_opt(kBegin1, kEnd1, GStart, GDest, system.maxGrad, system.maxSlew, segS.tSp, system.gradRasterTime);
                end
            else
                % initial gradient
                segS.tSp = GSspr.riseTime+GSspr.flatTime+GSspr.fallTime+seg.i1(k)*system.gradRasterTime;
                if((loop_m == 1)&&strcmp(T2prep,'on')), kBegin1 = [-segP.GXfac*spiral.kStart+dKA(k,1) -segP.GYfac*spiral.kStart];
                else
                    kBegin1 = [-spiral.kStart+dKA(k,1) 0];
                end
                kEnd1 = [kt(seg.ikseg(k,1),1) kt(seg.ikseg(k,1),2)+dKA(k,2)];
                GStart = [0 0];
                if(strcmp(spmode,'cont')||(k ~= 2*floor(k/2)))
                    GDest = [Gsp(seg.ikseg(k,1),1) Gsp(seg.ikseg(k,1),2)];
                end
                if(strcmp(spmode,'alt')&&(k == 2*floor(k/2)))
                    GDest = -[Gsp(seg.ikseg(k,1),1) Gsp(seg.ikseg(k,1),2)];
                end
                if(strcmp(segmode,'tan'))
                    [GBegin, tBegin, Ginit, tinit]  =  spiral_k2k_m(kBegin1, kEnd1, GStart, GDest, segS.tSp, system.maxSlew, system.gradRasterTime);
                else
                    if((loop_m == 1)&&strcmp(T2prep,'on')),[GBegin, tBegin,Gtran, tG, tt]  =  spiral_k2k_opt(kBegin1, kEnd1, GStart, GDest, system.maxGrad, system.maxSlew, segS.tSp-segP.t2, system.gradRasterTime);
                    else
                        [GBegin, tBegin,Gtran, tG, tt]  =  spiral_k2k_opt(kBegin1, kEnd1, GStart, GDest, system.maxGrad, system.maxSlew, segS.tSp, system.gradRasterTime);
                    end
                end
            end
            % terminal gradient
            ntk = acq.ntot-seg.i1(k)-seg.nkseg(k);
            segS.tSp = GSspr.riseTime+GSspr.flatTime+GSspr.fallTime+ntk*system.gradRasterTime;
            kBegin2 = [kt(seg.ikseg(k,2),1) kt(seg.ikseg(k,2),2)+dKE(k,2)];
            kEnd2 = [spiral.kStart+dKE(k,1) 0];
            GDest = [0 0];
            if(strcmp(spmode,'cont')||(k ~= 2*floor(k/2)))
                GStart = [Gsp(seg.ikseg(k,2),1) Gsp(seg.ikseg(k,2),2)];
            end
            if(strcmp(spmode,'alt')&&(k == 2*floor(k/2)))
                GStart = -[Gsp(seg.ikseg(k,2),1) Gsp(seg.ikseg(k,2),2)];
            end
            if(strcmp(segmode,'tan'))
                [GEnd, tEnd, Gfin, tfin]  =  spiral_k2k_m(kBegin2, kEnd2, GStart, GDest, segS.tSp, system.maxSlew, system.gradRasterTime);
            else
                [GEnd, tEnd, Gfin, tfin]  =  spiral_k2k_opt(kBegin2, kEnd2, GStart, GDest, system.maxGrad, system.maxSlew, segS.tSp, system.gradRasterTime);
            end
            % combine gradients
            if((k == 1)&&strcmp(initmode,'rev'))
                Gtot = [GBegin(1:end,:); -Gsp(nkseg1:-1:seg.ikseg(k,1),:); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),:); GEnd(1:end,:);[0 0]];
            else
                if((k == 1)&&strcmp(initmode,'outin'))
                    Gtot(:,1) = [GEnd(end:-1:1,1); Gsp(seg.ikseg(k,2):-1:seg.ikseg(k,1),1); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),1); GEnd(1:end,1);[0]];
                    Gtot(:,2) = [GEnd(end:-1:1,2); Gsp(seg.ikseg(k,2):-1:seg.ikseg(k,1),2); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),2); GEnd(1:end,2);[0]];
                else     
                    if((k == 1)&&strcmp(initmode,'sym'))
                        Gtot(:,1) = [GEnd(end:-1:1,1); Gsp(seg.ikseg(k,2):-1:seg.ikseg(k,1),1); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),1); GEnd(1:end,1);[0]];
                        Gtot(:,2) = [-GEnd(:,2); -Gsp(seg.ikseg(k,2):-1:seg.ikseg(k,1),2); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),2); GEnd(1:end,2);[0]];
                    else    
                        if(strcmp(spmode,'cont')||(k ~= 2*floor(k/2)))
                            Gtot = [GBegin(1:end,:); Gsp(seg.ikseg(k,1)+1:seg.ikseg(k,2),:); GEnd(1:end,:);[0 0]];
                        end
                        if(strcmp(spmode,'alt')&&(k == 2*floor(k/2)))
                            Gtot = [GBegin(1:end,:); -Gsp(seg.ikseg(k,1):-1:seg.ikseg(k,2)+1,:); GEnd(1:end,:);[0 0]];
                        end
                        ttot = [tBegin(1:end) [(1:seg.nkseg(k))*system.gradRasterTime] tEnd(2:end)];
                    end
                end
            end
            nGtot = length(Gtot);
            ttot = [1:nGtot]*system.gradRasterTime;
            % split Gtot
            if((loop_m == 1)&&strcmp(T2prep,'on'))
                nsp = round((GSspr.riseTime+GSspr.flatTime+GSspr.fallTime-segP.t2)/system.gradRasterTime);
            else
                nsp = round((GSspr.riseTime+GSspr.flatTime+GSspr.fallTime)/system.gradRasterTime);
            end
            Gsp1 = Gtot(1:nsp,:);
            Gspiral = Gtot(nsp+1:nsp+acq.ntot,:);
            Gsp2 = Gtot(nsp+acq.ntot+1:end,:);
            Gsp1_x = mr.makeArbitraryGrad('x','system',system,'waveform',Gsp1(:,1),'first',0,'last',mean(Gtot(nsp:nsp+1,1)));
            Gsp1_y = mr.makeArbitraryGrad('y','system',system,'waveform',Gsp1(:,2),'first',0,'last',mean(Gtot(nsp:nsp+1,2)));
            if(strcmp(segmode,'single'))
                for kadc = 1:1
                    Gspiral_x(kadc) = mr.makeArbitraryGrad('x','system',system,'maxSlew',4.2576e+11,'waveform',Gspiral(:,1));
                    Gspiral_y(kadc) = mr.makeArbitraryGrad('y','system',system,'maxSlew',4.2576e+11,'waveform',Gspiral(:,2));
                end
            else
                Gspiral_x = mr.makeArbitraryGrad('x','system',system,'maxSlew',4.2576e+11,'waveform',Gspiral(:,1),'first',Gsp1_x.last,'last',mean(Gtot(nsp+acq.ntot:nsp+acq.ntot+1,1)));
                Gspiral_y = mr.makeArbitraryGrad('y','system',system,'maxSlew',4.2576e+11,'waveform',Gspiral(:,2),'first',Gsp1_y.last,'last',mean(Gtot(nsp+acq.ntot:nsp+acq.ntot+1,2)));
            end
            Gsp2_x = mr.makeArbitraryGrad('x','system',system,'waveform',Gsp2(:,1),'first',Gspiral_x.last,'last',0);
            Gsp2_y = mr.makeArbitraryGrad('y','system',system,'waveform',Gsp2(:,2),'first',Gspiral_y.last,'last',0);
            rfref.signal = refenvelope*acqP.rflip(loop_m)/180;
            dOffset(indk) = 0;
            if(kconc == 2)
                dOffset(indk) = -acq.concB(indk).*dz^2;
            end
            rfref.phaseOffset = dOffset(indk) + rfref_phase-2*pi*rfref.freqOffset*mr.calcRfCenter(rfref);
            if (loop_m>1)
                seq.addBlock(GS4,rfref);
            else
                if(strcmp(T2prep,'off'))
                    seq.addBlock(GS4,rfref);
                end
            end
            if(strcmp(scanmode,'trim')&&(k == ntrim(loop_s)))
                seq.addBlock(GS5);
                seq.addBlock(GRtrim,adc);
                seq.addBlock(GS7);
            else
                if (Gsp1_x.shape_dur>GS5.shape_dur)
                    fprintf('\ntrajectory not compatible with sequence parameters. \nChange parameters and/or buy faster/stronger gradients\n');
                    return
                end
                Gsp1_x.first = 0;
                if(loop_m>1)
                    seq.addBlock(Gsp1_x,Gsp1_y,GS5);
                else
                    if(strcmp(T2prep,'off'))
                        seq.addBlock(Gsp1_x,Gsp1_y,GS5);
                    else
                        seq.addBlock(Gsp1_x,Gsp1_y);
                    end
                end
                if (loop_dummy>0)
                    if ((loop_m == 1)&& strcmp(T2prep,'on'))
                        seq.addBlock(Gspiral_x,Gspiral_y,GS5p2,adc);
                    else
                        seq.addBlock(Gspiral_x,Gspiral_y,adc);
                    end
                else
                    seq.addBlock(Gspiral_x,Gspiral_y);
                end
                if(Gsp2_x.shape_dur>GS7.shape_dur)
                    fprintf('\ntrajectory not compatible with sequence parameters. \nChange parameters and/or buy faster/stronger gradients\n');
                    return
                end
                seq.addBlock(Gsp2_x,Gsp2_y,GS7);
            end
        end
        seq.addBlock(GS4);
        seq.addBlock(GS5);
        seq.addBlock(delayTR);
    end
end

%% save all variables to SPITSE struct
SPITSE.acq    = acq;
SPITSE.acqP   = acqP;
SPITSE.seg    = seg;
SPITSE.segP   = segP;
SPITSE.spiral = spiral;
SPITSE.NEcho  = acqP.necho;

% save seq objects
obj_names = {'adc', 'delayTE1', 'delayTE2', 'delayTR', 'GRpre', 'GRpre_s', 'GRspoi_x', 'GRspoi_y', 'GRtrim', 'GS1', 'GS2', 'GS3', 'GS3p1', 'GS3p2', 'GS4', 'GS5', 'GS5p1', 'GS5p2', 'GS7', 'GS7p1', 'Gsp1_x', 'Gsp1_y', 'Gsp2_x', 'Gsp2_y', 'Gspiral_x', 'Gspiral_y', 'gz_fs', 'rf_fs', 'rfex', 'rfref'};
for j = 1:numel(obj_names)
    if exist(obj_names{j}, 'var')
        eval(['SPITSE.seq_objs.' obj_names{j} '=' obj_names{j} ';' ]);
    end
end

end
     

        
      