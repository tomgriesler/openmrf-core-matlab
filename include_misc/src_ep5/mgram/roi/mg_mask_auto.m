function mask = mg_mask_auto(image, otsu_factor)

if nargin<2
    otsu_factor = 1;
end

if ndims(image) == 3
    mask = zeros(size(image));
    for j=1:size(image,1)
        mask(j,:,:) = mg_mask_auto(squeeze(image(j,:,:)), otsu_factor);
    end
else

    % Convert complex image to magnitude
    image = abs(image);

    % Normalize the magnitude image to [0, 1] range
    image_norm = mat2gray(image);

    % Apply adaptive thresholding (Otsu's method)
    image_bin = image_norm > graythresh(image_norm)*otsu_factor;

    % Use edge detection for refining the mask
    edges = edge(image_bin, 'Canny');

    % Fill holes within the brain region (to avoid holes in the mask)
    filled_mask = imfill(edges, 'holes');

    % Remove small noise regions using area opening
    cleaned_mask = bwareaopen(filled_mask, round(0.01 * size(image,1)*size(image,2)) );

    % Fill any remaining holes and finalize the mask
    final_mask = imfill(cleaned_mask, 'holes');

    % Output the final binary mask
    mask = final_mask;

end
    
end

