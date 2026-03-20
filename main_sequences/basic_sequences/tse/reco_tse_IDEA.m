%% select file
clear
[study_name, study_path] = uigetfile('*.*', 'All Files (*.*)', pulseq_get_rawdata_path() );
study = fullfile(study_path, study_name);

twix_obj = mapVBVD(study, 'ignoreSeg', 'removeOS');
rawdata  = squeeze(twix_obj.image());
if ndims(rawdata==4)
    rawdata = squeeze(mean(rawdata,4));
end

rawdata           = permute(rawdata, [2,3,1]);
ImageArr_Coils    = kspace2image(rawdata);
[ImageArr, cmaps] = openadapt(ImageArr_Coils);

zero_params.onoff  = 1;
zero_params.radius = 6.0;
zero_params.factor = 2.0;

Images = mg_zero_filling(Images, zero_params);

xtv(Images)

clear twix_obj rawdata ImageArr_Coils;
