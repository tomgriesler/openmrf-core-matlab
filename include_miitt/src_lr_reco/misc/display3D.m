function imbig = display3D(myImage,rows,cols,clims,varargin)
% DISPLAY3D -- displays a 3D image in mosaic (tiled) view
% 
% Inputs
% 1. myImage -- (double 3D) the 3D image
% 2. rows -- (double 1x1) the number of rows in the mosiac view
% 3. cols -- (double 1x1) the number of columns in the mosiac view
% 4. optional arguments can be passed as 'Name',Value pairs
%       doNorm -- TRUE to normalize each 2D image between 0-1
%       fliplr -- TRUE to flip all images horizontally
%       flipud -- TRUE to flip all images vertically
% 
% Date: 3/16/2020
% Jesse Hamilton
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

if nargin < 2
    rows = ceil(sqrt(size(myImage,3)));
    cols = rows;
end

bFlipUD = false;
bFlipLR = false;
bReverse = false;
shiftz = 0;
donorm = true;
zoom = [];
if ~isempty(varargin)
    for i = 1:2:length(varargin)
        if strcmp(varargin{i},'reverse')
            bReverse = varargin{i+1};
        end
        if strcmp(varargin{i},'fliplr')
            bFlipLR = varargin{i+1};
        end
        if strcmp(varargin{i},'flipud')
            bFlipUD = varargin{i+1};
        end
        if strcmp(varargin{i},'shiftz')
            shiftz = varargin{i+1};
        end
        if strcmp(varargin{i},'doNorm')
            donorm = varargin{i+1};
        end
        if strcmp(varargin{i},'zoom');
            zoom = varargin{i+1};
        end
    end
end

[ny,nx,nz] = size(myImage);

if shiftz ~= 0
    myImage = circshift(myImage,[0 0 shiftz]);
end
if bReverse
    myImage = myImage(:,:,end:-1:1);
end

z = 0;
imbig = zeros(rows*ny,cols*nx,'single');
for i = 1:rows
    for j = 1:cols
        z = z+1;
        if z > size(myImage,3), continue; end
        im = squeeze(myImage(:,:,z));
        if bFlipLR
            im = fliplr(im);
        end
        if bFlipUD
            im = flipud(im);
        end
        if donorm
            im = mat2gray(im);
        end
        imbig((i-1)*ny+1:i*ny,(j-1)*nx+1:j*nx) = im;
    end
end

if nargin < 4
    clims = [min(imbig(:)) max(imbig(:))];
end
if isempty(zoom)
    imagesc(imbig,clims); axis image; axis off;
else
    imshow(imbig,clims,'initialmagnification',zoom); 
end