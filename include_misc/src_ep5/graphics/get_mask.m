function dataset = get_mask(varargin)
% function get_mask takes a bruker dataset as input and returns a bruker dataset as
% output containing a field dataset.mask
%
%
% FLAGS : 
%   dataset    : must be a valid MR dataset
%   mask_type  : 'noise_level' , 'contrast', 'intensity', 'roi', 'import'
%   intensity_level : a value between 0 and 1. Default to 0.3 
%                    Used with mask_type='intensity'
%                    Consider all voxels whose relative signal intensity > intensity_level
%   import         : set FLAG to 1 and pass existing mask as flag importedMAsk : i.e ....,'mask_type','import,'importedMask',myMask,...
%   force_mask : 1 yes, 0 no
%   slice : defaults to 1, the first slice of the dataset
%   automatic_mode  : 0 no, 1 yes : if 'automatic_mode' is ON, then the mask will be build 
%                     without user interaction if the mask type supports a non-interactive mode
warning('MR:M_FileCalled', 'MR:get_mask called')


%assign default values
DATASET = [];
mask_type = 'ROI';
location= pwd;
force_mask=1;
%reference_slice = 1;
save_data=-1;
intensity_level = 0.3;
slice = 1;
automatic_mode = 0;
parseVarargin;

% WARNING DATASET will get defined twice inside VARARGIN !!!
% so remove it before this happends
vararginExcludeList = {'dataset'};
reduceVarargin;
varargin = reducedVarargin;
% WARNING : should it be assigned without check
force_noise = 0;
% shortcut name
referenceSlice = squeeze(abs(dataset.data(slice,:,:)));


% check that the dataset is of good type
if ~strcmp(dataset.myReco.datasetType,'Slice_Phase_Read')
   warning('MR:NotImplemented', 'MR:get_mask. Only Slice_Phase_Read mode is implemented')
   return
end


if ~isfield(dataset,'mask')
   dataset.mask = [];
end

% check if existing mask should be overwriten/recalculated
if force_mask
   % reset mask value
   if isfield(dataset,'mask')
      dataset.mask = [];
   end
   % reset SNR reference value in case Mask is calculated based on an SNR value
   dataset.mask.snr_ref = -1;	
end
   
% we already checked if one forced a new mask.
% now we do nothing if mask exists and not empty
if isfield(dataset, 'mask') && isfield(dataset.mask, 'data') && ~isempty(dataset.mask.data)
   return
end

% now we start calculating a new mask

switch mask_type
   case 'import'
      % create mask based on voxel intensity
      dataset.mask.data = importedMask;
   case 'contrast'
      % create mask based on voxel intensity
      imagesc(abs(referenceSlice))
      dataset.mask.data = makemask;
   case 'intensity'
      % WARNING FIXME : here we always overwrite the existing dataset.mask.intensity_level
      dataset.mask.intensity_level = intensity_level;
      maxSlice = max(referenceSlice(:));
      if automatic_mode == 0
         % calculate mask with user interaction
         finished = false;
         while not(finished)
           prompt = {'Enter the relative minimum intensity  (0-1) to be part of mask' };
           name = ['Define minimum intensity: [0-1] ' ];
           numlines = 1;
           defaultanswer = {num2str(dataset.mask.intensity_level)};
           dataset.mask.intensity_level = str2num(char(inputdlg(prompt , name , numlines , defaultanswer)));
           % create mask: i.e. do not fit data set with SNR lower than limit
           ListOfMaskElements = find( referenceSlice > (dataset.mask.intensity_level * maxSlice));
           dataset.mask.data = zeros(size(referenceSlice));
           dataset.mask.data(ListOfMaskElements) = 1;
           imagesc(dataset.mask.data);
           ButtonName=questdlg('What is your wish?', ...
                         'Change the intensity level and recalculate mask', ...
                          'Recalculate mask','Continue','Continue');
            switch ButtonName
               case 'Recalculate mask'
                 finished = false;
                case 'Continue'
	                 finished = true;
            end
         end
      else
         % calculate mask without user interaction 
         ListOfMaskElements = find( referenceSlice > (intensity_level * maxSlice));
         dataset.mask.data = zeros(size(referenceSlice));
         dataset.mask.data(ListOfMaskElements) = 1;
      end
   case 'noise'
      % call noise. here
      myVarargin = {'force_noise',0,'+return_dataset' varargin{1:end}};
      dataset = get_noise('dataset',dataset, myVarargin{1:end} ) ;
      if dataset.mask.snr_ref<0
         % interactive mask construction
         finished = false;
         tmpMask = referenceSlice./dataset.noise.data ;
         while not(finished)
           prompt = {'Enter the minimum acceptable Signal to noise ratio: ' };
           name = ['Define SNR limit: Max(SNR) is ' ];
           numlines = 1;
           defaultanswer = {num2str(dataset.noise.snr*0.4)};
           dataset.mask.snr_ref = str2num(char(inputdlg(prompt , name , numlines , defaultanswer)));
           % create mask: i.e. do not fit data set with SNR lower than limit
           ListOfMaskElements = find( tmpMask > dataset.mask.snr_ref);
           dataset.mask.data = zeros(size(referenceSlice));
           dataset.mask.data(ListOfMaskElements) = 1;
           imagesc(dataset.mask.data);
           ButtonName=questdlg('What is your wish?', ...
                      'Change the SNR limit and recalculate mask', ...
                          'Recalculate mask','Continue','Continue');
            switch ButtonName
               case 'Recalculate mask'
                  finished = false;
                case 'Continue'
	              finished = true;
               end
         end
      else
            % calculate mask from existing noise level
            tmpMask =  referenceSlice./dataset.noise.data ;
            ListOfMaskElements = find( tmpMask > dataset.mask.snr_ref);
            dataset.mask.data = zeros(size(referenceSlice));
            dataset.mask.data(ListOfMaskElements) = 1;
      end
	  
   case 'roi'
      myVarargin = {varargin{1:end}};
      dataset = get_roi('dataset',dataset, myVarargin{1:end} ) ;
   case 'ROI'
      myVarargin = {varargin{1:end}};
      myTmp = get_roi('dataset',dataset, myVarargin{1:end} ) ;
      dataset.mask.data = myTmp.mask;

end
   