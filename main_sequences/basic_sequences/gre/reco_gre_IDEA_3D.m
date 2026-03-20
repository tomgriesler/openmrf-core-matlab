%% select file
clear
[study_name, study_path] = uigetfile('*.*', 'All Files (*.*)', pulseq_get_rawdata_path() );
study = fullfile(study_path, study_name);

twix_obj = mapVBVD(study, 'ignoreSeg', 'removeOS');
rawdata  = squeeze(twix_obj.image());

if ndims(rawdata)==5
    rawdata = squeeze(mean(rawdata,5));
end

rawdata  = permute(rawdata, [4, 2, 3, 1]);
rawdata1 = squeeze(rawdata(1,:,:,:));
rawdata2 = squeeze(rawdata(2,:,:,:));
rawdata3 = squeeze(rawdata(3,:,:,:));

temp_images = kspace2image(rawdata1);
temp_images = openadapt(temp_images);
temp_images = rot90(temp_images, -1);
Images(1,:,:) = temp_images(:,:);

temp_images = kspace2image(rawdata2);
temp_images = openadapt(temp_images);
temp_images = flipud(temp_images);
Images(2,:,:) = temp_images(:,:);

temp_images = kspace2image(rawdata3);
temp_images = openadapt(temp_images);
temp_images = rot90(temp_images, 1);
Images(3,:,:) = temp_images(:,:);

zero_params.onoff  = 1;
zero_params.radius = 6.0;
zero_params.factor = 2.0;

Images = mg_zero_filling(Images, zero_params);

xtv(Images)

clear temp_images rawdata1 rawdata2 rawdata3 twix_obj;
