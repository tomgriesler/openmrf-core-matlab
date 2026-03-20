function images = mg_zero_filling(images, zero_params)

% Version: Maximilian Gram, 21.03.2024
%          Maximilian Gram, 28.10.2025 -> new filter types

%% default zero filling parameters
if nargin<2
    zero_params = [];
end
if isempty(zero_params)
    zero_params.filter = 'kaiser';
end
if ~isfield(zero_params, 'filter')
    zero_params.filter = 'kaiser';
end
if ~isfield(zero_params, 'radius')
    if strcmp(zero_params.filter, 'kaiser')
        zero_params.radius = 6.0;
    end
    if strcmp(zero_params.filter, 'tukey')
        zero_params.radius = 0.5;
    end
end
if ~isfield(zero_params, 'factor')
    zero_params.factor = 2.0;
end

%% single image case
if ismatrix(images)
    temp = images;
    clear images;
    images(1,:,:) = temp;
    clear temp;
end

%% create filter window
[Nimages, ny, nx] = size(images);
Nx = round(nx * zero_params.factor);
Ny = round(ny * zero_params.factor);

switch zero_params.filter
    case 'kaiser'
        Wx = kaiser(Nx, zero_params.radius);
        Wy = kaiser(Ny, zero_params.radius);
    case 'hann'
        Wx = hann(Nx);
        Wy = hann(Ny);
    case 'tukey'
        Wx = tukeywin(Nx, zero_params.radius);
        Wy = tukeywin(Ny, zero_params.radius);
    case 'none'
        Wx = ones(Nx, 1);
        Wy = ones(Ny, 1);
end

W = Wy * Wx';
W = W / max(abs(W(:)));
W = permute(repmat(W, [1, 1, Nimages]), [3, 1, 2]);

%% zero filling
kspaces = zeros(Nimages, Ny, Nx, 'like', images);
ys = floor((Ny - ny)/2) + 1;  ye = ys + ny - 1;  % y start/end
xs = floor((Nx - nx)/2) + 1;  xe = xs + nx - 1;  % x start/end
kspaces(:, ys:ye, xs:xe) = image2kspace(images);

%% apply filter function
kspaces = kspaces .* W;

%% inverse fourier transform
images = kspace2image(kspaces);

%% single image case
images = squeeze(images);

end
