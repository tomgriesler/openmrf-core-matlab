function [ mask_fit, mask_fit_3D ] = mg_get_mask_fit( Image, mode_holes, n_3D)

% Version: Maximilian Gram, 21.03.2024

if ndims(Image)==3
    Image = squeeze(mean(abs(Image)));
end


if nargin<2
    mode_holes = 'small_holes';
end

if nargin<3
    n_3D = [];
end

figure(9673)
imagesc(real(Image))
axis image; axis off; colormap gray;
mask_fit = mg_makemask();

if strcmp(mode_holes, 'small_holes')
    filled     = imfill(mask_fit, 'holes');
    holes      = filled & ~mask_fit;
    bigholes   = bwareaopen(holes, 100);
    smallholes = holes & ~bigholes;
    mask_fit   = mask_fit | smallholes;
elseif strcmp(mode_holes, 'holes')
    mask_fit = imfill(mask_fit, 'holes');
end

imagesc(mask_fit);
axis image; axis off; colormap gray;

if ~isempty(n_3D)
    mask_fit_3D = permute(repmat(mask_fit, [1 1 n_3D]), [3 1 2]);
else
    mask_fit_3D = [];
end

end

