function [rawdata, f0] = PRESS_reco(study)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% study: enter path of study .dat
% study: or enter '0' for using dropdown select

%% import/sort kSpace rawdata
twix_obj = mapVBVD(study, 'ignoreSeg');
rawdata  = squeeze(twix_obj.image());
f0 = twix_obj.hdr.Meas.lFrequency;

end

