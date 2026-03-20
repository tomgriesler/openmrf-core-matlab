function [adc, adcTime, adcBW, adcNSamples, adcDwell] = SPI_calc_adc(adcTimeDesired, adcBWDesired, system)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

    %% calc adc
    adcNSamplesDesired = adcTimeDesired * adcBWDesired;
    adcNSamples        = ceil( adcNSamplesDesired / 64 ) * 64;
    adcDwell           = ceil( adcTimeDesired / adcNSamples / system.adcRasterTime ) * system.adcRasterTime;
    adcTime            = adcDwell * adcNSamples;
    adcBW              = 1 / adcDwell;
    adc                = mr.makeAdc(adcNSamples, 'Dwell', adcDwell, 'system', system);
    adc.phaseOffset    = 0;

    %% prevent adc sampling truncation 
    if adcNSamples > 2^14
        error('maximum number of adc samples exceeds 2^14 (16384) -> decrease bandwidth!');
    end
    
    %% disp params
    disp(' ');
    disp('----------------------- ADC ----------------------- ');
    disp([ 'N samples desired: '  num2str(adcNSamplesDesired)  '     N samples used: '           num2str(adcNSamples)                ]);
    disp([ 'adc time desired: '   num2str(adcTimeDesired*1e3,  '%.1f') 'ms     adc time used: '  num2str(adcTime*1e3,  '%.1f') 'ms'  ]);
    disp([ 'adc BW desired: '     num2str(adcBWDesired/1e3,    '%.1f') 'kHz    adc BW used: '    num2str(adcBW/1e3,    '%.1f') 'kHz' ]);
    disp('--------------------------------------------------- ');
    disp(' ');

end

