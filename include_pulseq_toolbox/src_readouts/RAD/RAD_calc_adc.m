function [adc, adcTime, adcNSamples, adcBW, adcDwell] = RAD_calc_adc(t_adc, n_samples, system)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% calc adc
adcTimeDesired     = t_adc;
adcNSamplesDesired = n_samples;
adcBWDesired       = 1 / (t_adc/n_samples);
adcTime            = system.gradRasterTime * round( adcTimeDesired / system.gradRasterTime ) - system.gradRasterTime;
adcNSamples        = ceil(adcNSamplesDesired/8)*8;
adcDwell           = round( adcTime / adcNSamples /100e-9 ) *100e-9;
adcBW              = 1/adcDwell;
adc                = mr.makeAdc(adcNSamples, 'Dwell', adcDwell); 
adc.delay          = system.adcDeadTime;
adc.phaseOffset    = 0; % changed during rf spoiling

%% disp params
disp(' ');
disp('----------------------- ADC ----------------------- ');
disp([ 'N samples desired: '  num2str(adcNSamplesDesired)  '     N samples used: '           num2str(adcNSamples)                ]);
disp([ 'adc time desired: '   num2str(adcTimeDesired*1e3,  '%.1f') 'ms     adc time used: '  num2str(adcTime*1e3,  '%.1f') 'ms'  ]);
disp([ 'adc BW desired: '     num2str(adcBWDesired/1e3,    '%.1f') 'kHz    adc BW used: '    num2str(adcBW/1e3,    '%.1f') 'kHz' ]);
disp('--------------------------------------------------- ');
disp(' ');

end

