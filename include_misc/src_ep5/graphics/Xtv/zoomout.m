function zoomout(handles,state)



if state == 1
    %axHandle=get(handles.figure1,'currentaxes');
    %set(axHandle,'Buttondownfcn','box')
    figpos=get(handles.figure1,'position');
    axpos=get(handles.axes1,'position');
    %set(handles.figure1,'Pointer','arrow')
    
    %cpt=round(get(axHandle,'Currentpoint'))
    set(handles.figure1,'WindowButtonMotionFcn', '','WindowButtondownFcn', '', ...
        'Interruptible', 'off', 'BusyAction', 'Queue');
    set(handles.axes1,'ButtondownFcn', '')
    set(handles.figure1,'currentaxes',handles.axes1)
    img=getimage(handles.axes1);
    ylim=get(handles.axes1,'ylim');
    xlim=get(handles.axes1,'xlim');
    
    clim = caxis;
    imHandle=tv(img,3,[],1);
    caxis(clim);
    %set(handles.axes1,'dataaspectratio',handles.aspect)
    %axis('on','square','tight')
    set(handles.axes1,'ylim',ylim,'xlim',xlim)
    [ny, nx]=size(img);
    y=ny/handles.datasize(1);
    x=nx/handles.datasize(2);
    xtick=[0:x-1]*handles.datasize(2);
    ytick=[0:y-1]*handles.datasize(1);
    set(gca,'xgrid','on','xcolor','b','tickdir','in','xtick',xtick,'xticklabel',[])
    set(gca,'ygrid','on','ycolor','b','tickdir','in','ytick',ytick,'yticklabel',[])
    set(gcf,'Pointer','custom','PointerShapeCData',handles.pointer,'PointerShapeHotSpot',[2 2])
    %set(gcf,'Pointer','botr')
    
    %axpos_new=axpos.*figpos
    
    
    
    %function box
    waitforbuttonpress, cpt=round(get(handles.axes1,'Currentpoint'));rectpos=rbbox;
    cptend=round(get(handles.axes1,'Currentpoint'));
    
    point1 = cpt(1,1:2);              % extract x and y
    point2 = cptend(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    offset = abs(point1-point2);         % and dimensions
    
    length = point2-point1;
    length(1,3) = 1;
 
    x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
    y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
    hold on
    axis manual
    plot(x,y,'r','linewidth',3)          % draw box around selected region
    hold off
    
    pause(1)   
    
    set(handles.axes1,'Xlim',[cpt(1,1) cptend(1,1)]);
    set(handles.axes1,'ylim',[cpt(1,2) cptend(1,2)]);
    
end

set(handles.axes1,'dataAspectRatio',[handles.FOV./handles.datasize]);
set(handles.figure1,'WindowButtonMotionFcn', 'impix(''PointerMotion'')')
sp_contrastimage(handles)


