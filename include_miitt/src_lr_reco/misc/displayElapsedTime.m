function [minutes,seconds]=displayElapsedTime(tstart,tend)
elapsed = etime(tend,tstart);
minutes = floor(elapsed/60);
seconds = elapsed - minutes*60;
fprintf('  Time elapsed ... %.0f min, %.1f sec\n',minutes,seconds)
end