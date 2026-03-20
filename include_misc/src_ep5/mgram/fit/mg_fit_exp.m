function [A, B, R2, RSS] = mg_fit_exp(xdata, ydata, A_start, B_start, mod_fit)

% Version: Maximilian Gram, 21.03.2024

%% fit modus
% 0: matlab toolbox
% 1: fminsearch(), no boundarie
% 2: fmincon()

if nargin<5
    mod_fit = 1;
end

% case for no start values
if ((nargin<4) || (nargin<3))
    A_start = max(ydata);
    B_start = xdata(end) * 0.75;
end

% use column vectors
if ( size(xdata,1)<size(xdata,2) )
    xdata = xdata';
end
if ( size(ydata,1)<size(ydata,2) )
    ydata = ydata';
end

%% matlab toolbox
if mod_fit == 0
    myexp    = fittype('a*exp(-x/b)');
    [f, gof] = fit(xdata, ydata, myexp, 'StartPoint', [A_start, B_start], 'lower', [0.1*A_start, 0.1*B_start], 'upper', [10*A_start, 10*B_start]);
    myval    = coeffvalues(f);
    A        = myval(1);
    B        = myval(2);
    RSS      = gof.sse;
    R2       = gof.rsquare;
end

%% fminsearch(), no boundaries
if mod_fit == 1
    opts     = optimset('Display','off');
    P_start  = [A_start, B_start];
    FITfun   = @(P) P(1) * exp( -xdata/P(2) );
    RSSfun   = @(P) sum((ydata - FITfun(P)).^2);
    P_fit    = fminsearch(RSSfun, P_start, opts);
    weights  = xdata / P_fit(2);
    weights  = (-tanh(weights-5) + 1)/2;
    RSSfun   = @(P) sum((ydata - FITfun(P)).^2 .* weights);
    P_fit    = fminsearch(RSSfun, P_start, opts);
    A        = P_fit(1);
    B        = P_fit(2);
    TSS      = sum( (ydata-mean(ydata)).^2 );
    RSS      = RSSfun(P_fit);
    R2       = 1 - RSS/TSS;        
end

%% fmincon()
if mod_fit == 2
    P_start = [A_start,    B_start];
    P_lower = [A_start/10, B_start/10];
    P_upper = [A_start*10, B_start*10];
    FITfun  = @(P) P(1) * exp( -xdata/P(2) );
    RSSfun  = @(P) sum((ydata - FITfun(P)).^2);
    options = optimoptions('fmincon','Display','off','Algorithm','sqp');
    P_fit   = fmincon( RSSfun, P_start, [], [], [], [], P_lower, P_upper, [], options );
    A       = P_fit(1);
    B       = P_fit(2);
    TSS     = sum((ydata - mean(ydata)).^2);
    RSS     = RSSfun(P_fit);
    R2      = 1 - RSS/TSS;
end

end

