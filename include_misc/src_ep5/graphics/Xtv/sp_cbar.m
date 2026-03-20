function hbar=sp_cbar(handles)
%axpos=get(gca,'position');
%kids=get(gcf,'children');

hbar=handles.axes5;
state=get(handles.axes5,'Visible');

if isequal(state,'on'),
 
%if isempty(hbar)==0,pos=get(hbar,'position');
%    delete(hbar),
%    hbar=axes('position',pos);
%else
%hbar=colorbar
%return
%end

t = caxis(handles.axes1);
d = (t(2) - t(1))/size(colormap(handles.axes1),1);
t = [t(1)+d/2  t(2)-d/2];

n = size(colormap(handles.axes1),1);

img = image([0 1],t,(1:n)','Parent',hbar);
set(hbar,...
        'Ydir','normal',...
        'YAxisLocation','right',...
        'xtick',[],...
        'tag','Colorbar')
set(gcf,'currentaxes',handles.axes1)
end
%else, return,end%set(haxes,'position',axpos)