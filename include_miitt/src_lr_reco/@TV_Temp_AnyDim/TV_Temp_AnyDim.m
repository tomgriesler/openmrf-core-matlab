function  res = TV_Temp_AnyDim(dim1)


if nargin < 1
    dim1 = 6;
end
% finite-differencing operator along contrast-enhancement dimension.
%
% Li Feng 2016

res.adjoint = 0;
res.dim1 = dim1;
res = class(res,'TV_Temp_AnyDim');

