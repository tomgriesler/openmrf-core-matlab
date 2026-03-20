function res = Waveletxyt(filterType, filterSize, wavScale, imageSize)
% res = Wavelet([Filtertype, [filterSize, [wavScale)
%
% implements a wavelet operator
%
% (c) Michael Lustig 2007

if nargin < 3
	wavScale = 1;
end

if nargin <2
	filterSize = 4;
end

if nargin < 1
	filterType = 'Daubechies';
end

res.adjoint = 0;
res.imageSize = imageSize;
res.qmf = MakeONFilter(filterType, filterSize);
res.wavScale = wavScale;
res = class(res,'Waveletxyt');
