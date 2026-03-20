%% mgram 12.01.24
% test espirit function for calculating coil sensitivities
clear

load brain_8ch
clearvars -except DATA;

images_coils = permute(ifft2c(DATA), [3,1,2]);
[cmaps, mask, im] = mg_espirit_cmaps(images_coils, 0.02, 0.95, 24, [6,6]);

xtv(cmaps)
figure(); imagesc(mask); axis image; colormap gray;
xtv(im)

% [im_open, cmaps_open] = openadapt(images_coils);
% xtv(im_open)
% xtv(cmaps_open)