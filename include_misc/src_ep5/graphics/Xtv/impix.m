function impix(arg1)

if ischar(arg1) % execute a callback
   switch lower(arg1)
   case 'pointermotion'
      PointerMotion;
   otherwise
      error(['Unknown IMPIXEL command or callback: ' arg1 '.']);
   end 

elseif ishandle(arg1)
   
   images = findobj(arg1, 'type', 'image');
 
   
    %guidata(arg1,handles)
    %global handles
    
   if isempty(images)
      return
   end
end



%figureHandle=arg1;
figureHandle=get(0,'CurrentFigure');
name=get(figureHandle,'Name');

if isequal(name,'xtv')==0, return,end
% weiss beim besten Willen nicht mehr was das soll!!!
%if isempty(findstr(num2str(figureHandle),'.')),return,end

images = findobj(figureHandle, 'type', 'image');


if isempty(images) 
   figureHandle=0; imageType=''; img=[]; x=0; y=0;
   return
end

set(figureHandle, 'Interruptible', 'off', 'busyaction', 'queue');

set(figureHandle, 'WindowButtonMotionFcn', 'impix(''PointerMotion'')',...
    'buttondownfcn',['sp_updatecontrast(' '''start''' ')'],...
    'Interruptible', 'on', 'BusyAction', 'Queue');
%set(figureHandle, 'ButtondownFcn','impix(''getpointer'')')
%set(figureHandle, 'KeyPressFcn', 'print -f figureHandle -v', ...
%    'Interruptible', 'off', 'BusyAction', 'Queue');

set(figureHandle,'HandleVisibility','on')
[imageHandle,imageType,img,x,y] = OverImage(figureHandle);

if isempty(img), set(figureHandle,'Pointer','arrow'),...
                  set(figureHandle,'HandleVisibility','callback'),...
                  set(figureHandle,'buttonDownfcn',''),return 
elseif strmatch(imageType,'indexed'),return 
        
else 

set(figureHandle,'Pointer','cross'),

set(figureHandle,'HandleVisibility','on')
%set(figureHandle, 'ButtondownFcn','impix(''PointerPress'')')


% if isempty(get(findobj(figureHandle,'tag','line1')));
%     h1 = line([0 0],[0 0])
%     set(h1,'tag','line1')
% end
%    
% 
% if isempty(get(findobj(figureHandle,'tag','line2')));
%     h2 = line([0 0],[0 0])
%     set(h2,'tag','line2')
% end
% 
% set(findobj(figureHandle,'tag','line1'),'XData',[x x],'YData',[1 size(img,1)]);
% set(findobj(figureHandle,'tag','line2'),'XData',[1 size(img,2)],'YData',[y y]);

handles=guidata(figureHandle);
handles.x=round(x);
handles.y=round(y);
handles.img=img;
guidata(figureHandle,handles)
plotprofiles(handles)
%set(handles.figure1,'HandleVisibility','Callback')
end
%set(handles.figure1,'HandleVisibility','off')

%-------------------------------------------------------------
% Subfunction PointerMotion
%-------------------------------------------------------------

function PointerMotion

imageHandle = gcbo;
figureHandle=get(0,'CurrentFigure');
figureHandle = get(get(figureHandle,'Parent'),'Parent');
set(figureHandle,'HandleVisibility','on')
images = findobj(figureHandle, 'type', 'image');

if isempty(images), return,end
%set(figureHandle, 'Pointer','cross')

function PointerPress
load zoompointer
figureHandle=get(0,'CurrentFigure');
set(figureHandle,'Pointer','custom','PointerShapeCData',zoompointer,'PointerShapeHotSpot',[16 16])
