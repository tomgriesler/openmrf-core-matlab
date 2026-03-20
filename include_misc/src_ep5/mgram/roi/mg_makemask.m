function [mask] = mg_makemask()

% Version: Maximilian Gram, 21.03.2024

% changes contrast in current figure using contrastimage
% output is mask
% 
% by fxbreuer 19.08.05

% fprintf('Please alter image contrast until contrast is black-white.\nPress RETURN to continue ... \n')
title('Please alter image contrast by dragging the mouse until contrast is black-white. Press RETURN to continue.')
set(gcf,'currentcharacter','d')
contrastimage(gcf);
%set(gcf,'selected','on')
figure(gcf)
waitfor(gcf,'currentCharacter',char(13))
figure(gcf)
lim=get(gca,'clim');
img=getimage(gca);
out=img;
out(out<lim(1))=0;
out(out>lim(1))=1;

mask=out;

end



