function [imageHandle,imageType,img,x,y] = OverImage(figHandle)
% Return the index of which image we are over, and return a 0 if we
% aren't above an image.
images = findobj(figHandle, 'type', 'image');
if isempty(images)
   imageHandle=0; imageType=''; img=[]; x=0; y=0;
   return
end
% Make sure that the Image's Button Down & Up functions will queue
%set(images, 'ButtonDownFcn', 'pixval(''ButtonDownOnImage'')', ...
%   'Interruptible', 'off', 'BusyAction', 'Queue');
axHandles = get(images, {'Parent'});
axPositions = get([axHandles{:}], {'Position'});
axCurPt = get([axHandles{:}], {'CurrentPoint'});

% Loop over the axes, see if we are above any of them
imageHandle = 0;  
for k=1:length(axHandles)
   XLim = get(axHandles{k}, 'XLim');
   YLim = get(axHandles{k}, 'YLim');
   pt = axCurPt{k};
   x = pt(1,1); y = pt(1,2);
   if x>=XLim(1) & x<=XLim(2) & y>=YLim(1) & y<=YLim(2)
      imageHandle = images(k);
      break;
   end
end
% Figure out image type
if imageHandle ~= 0
   [img,flag] = getimage(imageHandle);
   switch flag
   case 1
      imageType = 'indexed';
   case {2,3}  % Grayscale or binary
      imageType = 'intensity';
   case 4
      imageType = 'rgb';
   otherwise
      error(['Invalid image, GETIMAGE returned flag = ' flag '.']);
   end
else
   imageHandle=0; imageType=''; img=[]; x=0; y=0;
end

   