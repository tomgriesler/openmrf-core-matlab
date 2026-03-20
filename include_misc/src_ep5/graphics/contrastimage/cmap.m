function hbar=cmap(haxes)
axpos=get(gca,'position');
kids=get(gcf,'children');

hbar=findobj(kids,'tag','Colorbar');
if isempty(hbar)==0,pos=get(hbar,'position');
    delete(hbar),
    hbar=axes('position',pos);
else
%hbar=colorbar
return
end

t = caxis(haxes);
d = (t(2) - t(1))/size(colormap(haxes),1);
t = [t(1)+d/2  t(2)-d/2];

n = size(colormap(haxes),1);

img = image([0 1],t,(1:n)','Parent',hbar);
set(hbar,...
        'Ydir','normal',...
        'YAxisLocation','right',...
        'xtick',[],...
        'tag','Colorbar')
set(gcf,'currentaxes',haxes)
set(haxes,'position',axpos)