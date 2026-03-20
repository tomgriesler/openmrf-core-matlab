function [A, B, R2, CORR, A_error, B_error] = mg_fit_lin(xdata, ydata, A_start, B_start, mod_fit )

% Version: Maximilian Gram, 21.03.2024

%% fit modus
% 0: matlab toolbox
% 1: fminsearch(), no boundarie

if nargin<5
    mod_fit = 1;
end

% use column vectors
if ( size(xdata,1)<size(xdata,2) )
    xdata = xdata';
end
if ( size(ydata,1)<size(ydata,2) )
    ydata = ydata';
end

% case for no start values
if nargin<3
    B_start = (ydata(end)-ydata(1)) / (xdata(end)-xdata(1));
    A_start = mean(ydata) - mean(xdata)*B_start;
end

%% matlab toolbox
if mod_fit == 0
    mylin    = fittype('a+b*x');
    [f, gof] = fit(xdata, ydata, mylin, 'StartPoint', [A_start, B_start], 'lower', [0.1*A_start, 0.1*B_start], 'upper', [10*A_start, 10*B_start]);
    myval    = coeffvalues(f);
    A        = myval(1);
    B        = myval(2);
    R2       = gof.rsquare;
    ci       = confint(f, 0.6827);
    A_error  = abs(ci(1,1)-ci(2,1))/2;
    B_error  = abs(ci(1,2)-ci(2,2))/2;
    CORR     = mg_get_pearson_correlation(xdata, ydata);    
end

%% fminsearch(), no boundaries
if mod_fit == 1
    P_start = [A_start, B_start];
    FITfun  = @(P) P(1) + P(2).*xdata;
    RSSfun  = @(P) sum((ydata - FITfun(P)).^2);            
    P_fit   = fminsearch(RSSfun, P_start);
    A       = P_fit(1);
    B       = P_fit(2);
    TSS     = sum( (ydata-mean(ydata)).^2 );
    RSS     = RSSfun(P_fit);
    R2      = 1 - RSS/TSS;
    CORR    = mg_get_pearson_correlation(xdata, ydata); 
end

end

