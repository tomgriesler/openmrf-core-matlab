function [rawdata, NImages, NRead, NCoils, rawdata_noise] = SPI_get_rawdata(twix_obj, Nnoise)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

if nargin<2
    Nnoise = 0;
end

%% import/sort kSpace rawdata
rawdata10D = twix_obj.image.unsorted();
nSegments  = twix_obj.image.NSet;

% force 10D format using dummy dimensions...
if ndims(rawdata10D)<3 % for the case nImages=1 and nSegments=1 
    temp_raw = rawdata10D;
	rawdata10D = zeros( size(temp_raw,1), size(temp_raw,2), 1, 1, 1, 1, 1, 1, 1, 2);
    rawdata10D(:,:,1,1,1,1,1,1,1,1) = temp_raw(:,:);
end
if ndims(rawdata10D)<4 % for the case nImages>1 and nSegments=1
    temp_raw = rawdata10D;
	rawdata10D = zeros( size(temp_raw,1), size(temp_raw,2), size(temp_raw,3), 1, 1, 1, 1, 1, 1, 2);
    rawdata10D(:,:,:,1,1,1,1,1,1,1) = temp_raw(:,:,:);
end

NRead   = size(rawdata10D,1);
NCoils  = size(rawdata10D,2);
NImages = size(rawdata10D,3);
rawdata = zeros( NImages, NRead*nSegments, NCoils );

if nSegments>1 
    for j=1:NCoils
    for k=1:NImages
    for l=1:nSegments
        temp(:) = rawdata10D(:,j,k,1,1,1,1,1,1,l);
        rawdata( k, (l-1)*NRead+1:l*NRead, j ) = temp(:);        
    end 
    end
    end
    NRead = NRead*nSegments;    
else
    for j=1:NCoils
    for k=1:NImages
        temp(:) = rawdata10D(:,j,k,1,1,1,1,1,1,1);
        rawdata( k, :, j ) = temp(:);        
    end
    end
end

if Nnoise>0
    rawdata_noise = rawdata(1:Nnoise,:,:);
    rawdata(1:Nnoise,:,:) = [];
    NImages = NImages - Nnoise;
else
    rawdata_noise = [];
end

end