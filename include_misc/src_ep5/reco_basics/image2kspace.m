function kspace = image2kspace(myImage, varargin)
% image2kspace    can perform 1D or 2D transformation from image to kspace 
%
% kspace = image2kspace(myImage,vector)
%      when vector is a vector of elements with value 0 or 1 with as 
%      many elements as the number of dimensions of myImage, then the transformations
%      are applied along dimensions when 1.
% kspace = image2kspace(myImage)
%     KSPACE can be a 2D matrix or a Multidimensional array of structure [images ... images Phase Read]
%     The transformation is a 2D transformation along the PHASE and READ dimensions only
%
% See image2kspace

%warning('MR:M_FileIterCalled', 'MR::image2kspace called')

% parse varargin input
v = [];
if nargin>1
    v = varargin{1};
end

dim = size(myImage);
% check if user indicated which dimensions to transform
if isempty(v)
   % assumes [images ... images Phase Read] structure
   v = zeros(1,length(dim));
   v(end-1:end) = 1; 
elseif length(v)<length(dim)
    % check special case were kspace is a 1D vector
    if isvector(myImage) && length(v)==1
        v = [v v];
    elseif length(v)==1 && v==1
        % assumes [images ... images Phase Read] structure and fft along Read encoding
        v = zeros(1,length(dim));
        v(end) = 1;
    elseif length(v)==1 && v==2
        % assumes [images ... images Phase Read] structure and fft along first phase encoding
        v = zeros(1,length(dim));
        v(end-1) = 1;
    else
        warning('MR::kspace2image bad size. Expending ')
        v(length(dim))=0;
    end
end

 
kspace = myImage;
for k=1:numel(dim)
    if v(k)==1
        % check if FFT necessary
        if dim(k) > 1
       %  zero-frequency component is assumed to be in center of spectrum.
            % but fft or ifft assume that the zero-frequency components is
            % at begining of spectrum
            % therefore we must first bring zero-frequency component to
            % beginning of "spectrum" for wanted dimension only
            % Because fftshift does not work for odd dimensions, one requires the use of ifftshift
            % which works for even and odds vectors like [-3 -2 -1 0 1 2 3]
            % or [-4 -3 -2 -1 0 1 2 3]
            % After Fourier transform, one would like to bring the
            % zero-frequency/time/whatever back to center of "spectrum"     
            kspace = fftshift(fft(ifftshift(kspace,k),[],k),k);
       end
    end
end
   
