function contrastimage(figHandle)

%
%   Function to provide dynamic adjustment of image contrast
%
%   Once the image is displayed, you can hold down a mouse button
%   and move the mouse to change the image contrast.
%
%   Needs helper function updatecontrast.m
%
%   Mark Griswold 

if nargin==1,figHandle=figHandle;else figHandle=gcf;end


images = findobj(figHandle, 'type', 'image');
if isempty(images) ,return,end
im=getimage(images(1));
if isa(im,'uint8'),return,end
[ny,nx]=size(im);
contstruct.clim=[min(im(:)) max(im(:))];

if isempty(get(findobj(gca,'tag','contdata')))
	rmax0=max(max(im));
	rmin0=min(min(im));
    %rmin0=0;
	contstruct.windowval=0.5;%*(rmax0-rmin0);
	contstruct.centerval=contstruct.windowval/2;
    contstruct.cx=0.5;
    contstruct.wx=0.5;
	contstruct.startpoint=1;
    contstruct.lastpoint=[ny/2 nx/2 1; ny/2 nx/2 0];
	contstruct.scaley=1./(5*ny);
    contstruct.scalex=1./(5*nx);%(rmax0-rmin0)./100;
	axtag=get(gca,'tag');
    contsruct.ax=gca;
	%cla
	%imagesc(im);
	%axis('square');
    set(gca,'tag',axtag);
	contstruct.axis=get(gca,'children');
    %contstruct.clim=get(gca,'clim');
 
	t1=text(10,10,' ','tag','contdata','userdata',contstruct,'color','w','fontsize',8);
    %bdf=get(image,'buttondownfcn')
	set(contstruct.axis,'buttondownfcn',['updatecontrast(' '''start''' ')']);
	%set(gca,'xgrid','off','xtick',[],'ytick',[]);
    set(gcf,'doublebuffer','on');
else
	axtag=get(gca,'tag');
	contstruct=get(findobj(gca,'tag','contdata'),'userdata');
	[n,m]=size(contstruct.axis);
	[yy,xx]=size(im);
    %set(contstruct.axis,'buttondownfcn',['updateimage3d(' '''start''' ')']);
	set(contstruct.axis(n),'cdata',im);
	kids=get(gca,'children');
	set(kids,'visible','off');
	set(contstruct.axis(n),'visible','on');
	set(gca,'ylim',[0.5 yy+0.5],'xlim',[0.5 xx+0.5]);
        set(gca,'tag',axtag);
%	set(findobj(gca,'tag','contdata'),'userdata',contstruct);
	set(gca,'xtick',[],'ytick',[]);
end






