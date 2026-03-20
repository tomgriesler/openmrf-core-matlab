function [bw,x,y] = roi( varargin )
%  ROI    : defines Region Of Interest
%  from Felix Breuer, imported 2005/08/02

% parameters
imageHandle  = [];
axisHandle   = [];
figureHandle = [];
img          = [];

if ( isnumeric(varargin{1}) )
   img = varargin{1};
   varargin = varargin(2:end);
end

% classic definition of varargin parameters
parseVarargin;

% make use of correct figures and axes
if ~(isempty(figureHandle))
   figure(figureHandle);
end

if isempty(imageHandle) 
   if isempty(axisHandle)
      imageHandle = imagesc(abs(img));
   else
      subplot(axisHandle)
      imageHandle = imagesc(abs(img));
   end
else
    if isempty(axisHandle)
       axisHandle = get(imageHandle, 'Parent');
      subplot(axisHandle)
    else
       subplot(axisHandle) % whatever imageHandle parent is
    end
end


%  if ~isempty(imageHandle) 
%     axisHandle = get(imageHandle, 'Parent');
%     figureHandle = get(axisHandle, 'Parent');
%  else
%     if ~isempty(axisHandle)
%        figureHandle = get(axisHandle, 'Parent');
%        figure(figureHandle)
%        subplot(axisHandle)
%        imageHandle = imagesc(abs(img)); % img must exist !!
%      else
%      
%      end
%     end
%  
%  end



[x,y]=getpoints( imageHandle );
bw=roipoly(img,x,y);





function varargout=getpoints(figHandle)

    
    if nargin==1,figHandle=figHandle;else figHandle=gcf;end


    images = findobj(figHandle, 'type', 'image');
    if isempty(images) ,return,end
    
  if isempty(get(gca,'userdata')),   
     contstruct.x=[];
     contstruct.y=[];
     axtag=get(gca,'tag');
     contsruct.ax=gca;
     contstruct.ima=getimage(images(1));
     %set(gca,'tag',axtag);
     contstruct.axis=get(gca,'children');
     set(gca,'userdata',contstruct);
     set(gca,'tag','uncompleted')
     set(contstruct.axis,'buttondownfcn',['updatepoints(' '''start''' ')']);
     set(gcf,'doublebuffer','on');
     %get(contstruct.axis,'ButtondownFcn')
     %waitfor(gcf,'selectiontype','open')    
     waitfor(contstruct.axis,'ButtondownFcn','')
     %set(gcf,'selectiontype','normal')
    % Return the answer
     contstruct=get(gca,'userdata');
     %set(contstruct.axis,'buttondownfcn','');
     x=contstruct.x(1:end);
     y=contstruct.y(1:end);
     set(gca,'userdata','')
    if (nargout >= 2)
        varargout{1} = x;
        varargout{2} = y;
    else
        % Grandfathered output syntax
        varargout{1} = [x(:) y(:)];
    end

 else    
       set(gca,'userdata','')
 end
 
 
 
 
% helperfunction 
 
function updatepoints(action)
% global windowval centerval startpoint startx starty imsize currentgamma cx wx
switch action

case 'start'
    mouse=get(gcf,'selectiontype');
    if strcmp('normal',mouse)    
    contstruct=get(gca,'userdata');
   	set(gcbf,'WindowButtonMotionFcn','updatepoints move')
   	set(gcbf,'WindowButtonUpFcn','updatepoints stop')
	
    set(gcbf,'pointer','crosshair')
    end
  if strcmp('alt',mouse)
    updatepoints close
  end    
   if strcmp('extend',mouse)
   
   %set(gcbf,'WindowButtondownFcn','updatepoints remove')  
   updatepoints remove
  end
    %hbar=colorbar;set(hbar,'ylim',contstruct.clim);    
    %end
    
    updatepoints move

    
case 'move'
    set(gcbf,'pointer','crosshair')
    contstruct=get(gca,'userdata');
    x=contstruct.x;
    y=contstruct.y;
    points=get(gca,'Currentpoint');
    xi=(points(1,1));
    yi=(points(1,2));
    if size(x,1)>1
    h=line(x(end-1:end),y(end-1:end));
    set(h,'color','yellow')
    end
    %[yi,xi]=[points(1,1),points(1,2)]
    x=[x;xi];
    y=[y;yi];
    contstruct.x=x;
    contstruct.y=y;
    
    set(gca,'userdata',contstruct)
    
case 'stop'
    contstruct=get(gca,'userdata');
    set(gcf,'windowbuttondownfcn','updatepoints start')
    set(gcf,'keypressFcn','updatepoints remove')
    set(gcf,'WindowButtonMotionFcn','')
   	set(gcf,'WindowButtonUpFcn','')
    %set(gcf,'windowbuttondownfcn','updatepoints start')
case 'close'
     contstruct=get(gca,'userdata');
     x=contstruct.x(1:end-1);
     y=contstruct.y(1:end-1);
    h=line([x(end);x(1)],[y(end);y(1)]);
    set(h,'color','yellow')
    set(gcf,'WindowButtonMotionFcn','')
    set(gcf,'WindowButtondownFcn','')
    contstruct.x=x;
    contstruct.y=y;
    %get(contstruct.axis,'ButtondownFcn')
    set(gca,'userdata',contstruct)
    set(contstruct.axis,'ButtondownFcn','')
    
        set(gcf,'pointer','arrow')
    %ima=contstruct.ima;
    %out=roipoly(ima,x,y);
    %imagesc(out)
    
case 'remove'
    set(gcf,'pointer','circle')
    contstruct=get(gca,'userdata');
    key=get(gcf,'currentcharacter');
    if or(isequal(key,char(8)), isequal(key,char(127))),
        x=contstruct.x(1:end-1); y=contstruct.y(1:end-1); 
        contstruct.x=x;contstruct.y=y; 
        
       lines=findobj(gca,'type','line'); delete(lines(1)); 
        set(gca,'userdata',contstruct)
    end
end

