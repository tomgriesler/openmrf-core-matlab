function [fit_results] = mg_fit_sin_complex(ImageArr, mask_fit, phase_arr)

% Version: Maximilian Gram, 21.03.2024

%% init
Np = size(ImageArr,1); % number of phases
Ny = size(ImageArr,2); % Nxy
Nx = size(ImageArr,3); % Nxy

Offset_Map     = zeros(Nx,Nx);
Amplitude_Map  = zeros(Ny,Nx);
Phase_Map      = zeros(Ny,Nx);
Modulation_Map = zeros(Ny,Nx);
R2_Map         = zeros(Ny,Nx);

Amplitude_STD_Map = zeros(Ny,Nx);

%% create fit mask (if not existent)
if nargin<2
    mask_fit = mg_get_mask_fit(squeeze(mean(abs(ImageArr))));
end 

%% sin fit
if nargin<3
    xdata = linspace(0, 2*pi, Np+1)';
    xdata(end) = [];
else
    xdata = phase_arr;
end

parfor y=1:Ny
for    x=1:Nx
if     mask_fit(y,x) == 1

    % fit data
    ydata    = squeeze(ImageArr(:,y,x));
    ydata_re = real(ydata);
    ydata_im = imag(ydata);
    
    % fit model
    % y_c  = y_re + 1i * y_im
    % y_re = P(1) + P(3) * sin( P(5) + P(6) * x )
    % y_im = P(2) + P(4) * sin( P(5) + P(6) * x )
    
    % start params for 1st fit
    P1_start = mean(ydata_re);
    P2_start = mean(ydata_im);    
    P3_start = (max(ydata_re)-min(ydata_re)) / 2;
    P4_start = (max(ydata_im)-min(ydata_im)) / 2;
    P5_start = 0;
    
    %
    options = optimset('Display','off');

    % fminsearch: 1st fit          
    P_start   = [ P1_start, P2_start, P3_start, P4_start, P5_start ];
    FITfun_re = @(P) P(1) + P(3) * sin( P(5) + 1 * xdata );
    FITfun_im = @(P) P(2) + P(4) * sin( P(5) + 1 * xdata );
    RSSfun_re = @(P) sum( (ydata_re - FITfun_re(P)).^2 );
    RSSfun_im = @(P) sum( (ydata_im - FITfun_im(P)).^2 );
    RSSfun    = @(P) RSSfun_re(P) + RSSfun_im(P);
    P_fit     = fminsearch(RSSfun, P_start, options);
    P_fit(5)  = wrapTo2Pi(P_fit(5));
    
    % fminsearch: 2nd fit
    P1_start  = P_fit(1);
    P2_start  = P_fit(2);
    P3_start  = P_fit(3);    
    P4_start  = P_fit(4);
    P5_start  = P_fit(5);
    P6_start  = 1;
    P_start   = [ P1_start, P2_start, P3_start, P4_start, P5_start, P6_start ];
    FITfun_re = @(P) P(1) + P(3) * sin( P(5) + P(6) * xdata );
    FITfun_im = @(P) P(2) + P(4) * sin( P(5) + P(6) * xdata );
    RSSfun_re = @(P) sum( (ydata_re - FITfun_re(P)).^2 );
    RSSfun_im = @(P) sum( (ydata_im - FITfun_im(P)).^2 );
    RSSfun    = @(P) RSSfun_re(P) + RSSfun_im(P);    
    P_fit     = fminsearch(RSSfun, P_start, options);
    P_fit(5)  = wrapTo2Pi(P_fit(5));

    % fminsearch: 3rd fit
    P1_start  = P_fit(1);
    P2_start  = P_fit(2);
    P3_start  = P_fit(3);    
    P4_start  = P_fit(4);
    P5_start  = P_fit(5);
    P6_start  = P_fit(6);
    P_start   = [ P1_start, P2_start, P3_start, P4_start, P5_start, P6_start ];
    FITfun_re = @(P) P(1) + P(3) * sin( P(5) + P(6) * xdata );
    FITfun_im = @(P) P(2) + P(4) * sin( P(5) + P(6) * xdata );
    RSSfun_re = @(P) sum( (ydata_re - FITfun_re(P)).^2 );
    RSSfun_im = @(P) sum( (ydata_im - FITfun_im(P)).^2 );
    RSSfun    = @(P) RSSfun_re(P) + RSSfun_im(P);    
    P_fit     = fminsearch(RSSfun, P_start, options);
    P_fit(5)  = wrapTo2Pi(P_fit(5));
    
    % invert phase and amplitude
    if P_fit(5)>pi
        P_fit(3) = -P_fit(3);
        P_fit(4) = -P_fit(4);
        P_fit(5) = P_fit(5) - pi;
    end
    
    Offset_Map(y,x)     = P_fit(1) + 1i * P_fit(2);
    Amplitude_Map(y,x)  = P_fit(3) + 1i * P_fit(4);
    Phase_Map(y,x)      = P_fit(5);
    Modulation_Map(y,x) = P_fit(6);
    
    % calculate R2
    TSS_re      = sum( (ydata_re - mean(ydata_re)).^2 );
    TSS_im      = sum( (ydata_im - mean(ydata_im)).^2 );
    TSS         = TSS_re + TSS_im;       
    R2_Map(y,x) = 1 - RSSfun(P_fit) / TSS;
    
    % calculate amplitde map based on std()
    Amplitude_STD_Map(y,x) = std(ydata);

end
end
end

%% output data
fit_results.Offset_Map        = Offset_Map;
fit_results.Amplitude_Map     = Amplitude_Map;
fit_results.Phase_Map         = Phase_Map;
fit_results.Modulation_Map    = Modulation_Map;
fit_results.R2_Map            = R2_Map;
fit_results.Amplitude_STD_Map = Amplitude_STD_Map;

end

