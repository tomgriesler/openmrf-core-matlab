%% select file
clear
[study_name, study_path] = uigetfile('*.*', 'All Files (*.*)' );
study = fullfile(study_path, study_name);

twix_obj = mapVBVD(study, 'ignoreSeg', 'removeOS');
twix_obj = twix_obj{2};
rawdata  = squeeze(twix_obj.image());

if ndims(rawdata)==4
    rawdata = squeeze(mean(rawdata,4));
end

rawdata      = permute(rawdata, [2, 1, 3]);
Images_coils = kspace2image(rawdata);
Images       = openadapt(Images_coils);

zero_params.onoff  = 1;
zero_params.radius = 6.0;
zero_params.factor = 2.0;

Images = mg_zero_filling(Images, zero_params);

xtv(Images)

clear Images_coils rawdata twix_obj;
