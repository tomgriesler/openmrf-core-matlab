function plotprofiles(Handle)


%if size(Handle.img,1)>Handle.datasize(1),return, end
%if size(Handle.img,2)>Handle.datasize(2),return, end


ny=round(size(Handle.img,1)/Handle.datasize(1));
nx=round(size(Handle.img,2)/Handle.datasize(2));

imgnrl=ceil(Handle.x/Handle.datasize(2));
imgnrr=ceil(Handle.y/Handle.datasize(1));
imgnr=(imgnrr-1)*nx+imgnrl;
x=mod(Handle.x,Handle.datasize(2));
y=mod(Handle.y,Handle.datasize(1));
if imgnr<=Handle.datasize(3), z=imgnr; else z=0; return, end
if y==0, y=Handle.datasize(1);end
if x==0, x=Handle.datasize(2);end
if size(Handle.img,2)==Handle.datasize(2) & ... 
   size(Handle.img,1)==Handle.datasize(1), z=Handle.z;
end

% 

% 
axes(Handle.axes1)

updatestr = sprintf('x= %g,   y= %g,   NR.= %g,  ==> %6.4f',x,y,z,Handle.data(y,x,z));
set(Handle.update,'String',updatestr)
set(Handle.figure1,'CurrentAxes',Handle.axes2)
xlim=get(Handle.axes1,'ylim');
if xlim(2)>Handle.datasize(1),xlim=[0 Handle.datasize(1)];end
set(Handle.axes2,'xlim',xlim)
plot(squeeze(Handle.data(:,x,z)),'k')
stepx=5*round(round(Handle.datasize(1)/10)./5);
set(Handle.axes2,'Xlim',[0 Handle.datasize(1)],'Xtick',[0:stepx:Handle.datasize(1)])
        
stdval=std(squeeze(Handle.data(:,x,z)));
maxval=max(squeeze(Handle.data(:,x,z)));
minval=min(squeeze(Handle.data(:,x,z)));



if isequal(minval,maxval) | isnan(minval) | isnan(maxval) | isinf(minval) | isinf(maxval), 
set(Handle.axes2,'Ylimmode','auto');
else
set(Handle.axes2,'Ylim',[minval-2*stdval maxval+2*stdval]);
end
ylim=get(Handle.axes2,'Ylim');
line([y y],[ylim(1) ylim(2)])
set(Handle.text9,'string','U->D profile')



set(Handle.figure1,'CurrentAxes',Handle.axes3)
xlim=get(Handle.axes1,'xlim');
if xlim(2)>Handle.datasize(2),xlim=[0 Handle.datasize(2)];end
%set(Handle.axes3,'Xlim',xlim)
plot(squeeze(Handle.data(y,:,z)),'k')
set(Handle.axes3,'Ylimmode','auto');
stepx=5*round(round(Handle.datasize(2)/10)./5);
set(Handle.axes3,'Xlim',[0 Handle.datasize(2)],'Xtick',[0:stepx:Handle.datasize(2)])
        
stdval=std(squeeze(Handle.data(y,:,z)));
maxval=max(squeeze(Handle.data(y,:,z)));
minval=min(squeeze(Handle.data(y,:,z)));

if isequal(minval,maxval) | isnan(minval) | isnan(maxval) | isinf(minval) | isinf(maxval),
set(Handle.axes3,'Ylimmode','auto');
else
set(Handle.axes3,'Ylim',[minval-2*stdval maxval+2*stdval]);
end
ylim=get(Handle.axes3,'Ylim');
line([x x],[ylim(1) ylim(2)])
set(Handle.text10,'string','L->R profile')
set(Handle.figure1,'CurrentAxes',Handle.axes4)

if size(Handle.data,3)==1
    hist(Handle.data(find(Handle.data~=0)),100)
    set(Handle.text11,'string','histogram')
else
    stepx=2*round(round(Handle.datasize(3)/10)/2);
    %set(Handle.axes4,'Xlim',[1 Handle.datasize(3)],'Xtick',[1:stepx:Handle.datasize(3)])
    maxval=max(squeeze(Handle.data(y,x,:)));
    minval=min(squeeze(Handle.data(y,x,:)));
    if strcmp(num2str(maxval),num2str(minval)),
        set(Handle.axes4,'Ylim',[maxval-1 maxval+1]);
        plot(squeeze(Handle.data(y,x,:)),'k');
        %ylim=get(Handle.axes4,'Ylim')
    set(Handle.axes4,'Xlim',[0 Handle.datasize(3)],'Xtick',[0:stepx:Handle.datasize(3)])    
    else
        
        set(Handle.axes4,'Ylimmode','auto');
        
        stdval=std(squeeze(Handle.data(y,x,:)));
        maxval=max(squeeze(Handle.data(y,x,:)));
        minval=min(squeeze(Handle.data(y,x,:)));
        set(Handle.axes4,'Ylim',[minval-2*stdval maxval+2*stdval]);
        plot(squeeze(Handle.data(y,x,:)),'kx:')
        set(Handle.text11,'string','F->B profile')
        %set(Handle.axes4,'Ylim',[min(squeeze(Handle.data(y,x,1,:))) max(squeeze(Handle.data(y,x,1,:)))]);
        set(Handle.axes4,'Xlim',[0 Handle.datasize(3)],'Xtick',[0:stepx:Handle.datasize(3)])  
    end
    ylim=get(Handle.axes4,'Ylim');
   line([z z],[ylim(1) ylim(2)])
   
end
    %axes(Handle.axes1)
    
    
    %set(Handle.axes1,'HandleVisibility','Callback')
    %set(Handle.axes2,'HandleVisibility','off')
    %set(Handle.axes3,'HandleVisibility','off')
    %set(Handle.axes4,'HandleVisibility','off')
    
