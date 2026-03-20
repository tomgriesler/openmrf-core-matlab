% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% add dummy MRF preparations

% init loop counters for contrast preparations
loop_SAT    = 1;
loop_INV    = 1;
loop_T2     = 1;
loop_SL     = 1;
loop_MLEV   = 1;
loop_ADIASL = 1;
  
for loop_MRF = 1 : MRF.n_segm

    % reset seq struct
    clear seq;
    seq = mr.Sequence(system);
 
    % MRF Preparations:
    MRF_add_preparation();

    % calc prep durations
    MRF.prep_durations(loop_MRF,1) = seq.duration;

end
clear seq;
seq = mr.Sequence(system);

%% calculate dynamic delays
MRF.prep_max = max(MRF.prep_durations);
for loop_MRF = 1 : MRF.n_segm
    MRF.delay_dynamic(loop_MRF,1) = mr.makeDelay( round((MRF.prep_max-MRF.prep_durations(loop_MRF))/system.gradRasterTime)*system.gradRasterTime + system.gradRasterTime );
end