function updatecontrast(action)
%global windowval centerval startpoint startx starty imsize currentgamma cx wx
switch action

case 'start'
    mouse=get(gcf,'selectiontype');
    if strcmp('normal',mouse)    
   	set(gcbf,'WindowButtonMotionFcn','updatecontrast move')
   	set(gcbf,'WindowButtonUpFcn','updatecontrast stop')
	set(gcbf,'pointer','fleur')
    end
    if strcmp('alt',mouse)
    updatecontrast reset
    end    
    if strcmp('extend',mouse)
    set(gcbf,'WindowButtonUpFcn','')
    %hbar=colorbar;set(hbar,'ylim',contstruct.clim);    
    end
    
    updatecontrast move

case 'move'
	contstruct=get(findobj(gca,'tag','contdata'),'userdata');
   
         
	if contstruct.startpoint == 1
	  
		contstruct.firstpoint = get(gca,'currentpoint');
		contstruct.startpoint= 0;
		contstruct.startx=contstruct.firstpoint(1,1);
		contstruct.starty=contstruct.firstpoint(1,2);
		set(findobj(gca,'tag','contdata'),'userdata',contstruct);
		%set(findobj(gca,'tag','contdata'),'string',sprintf('C: %d W: %d',round(contstruct.centerval),round(contstruct.windowval)));
        %elseif contstruct.startpoint==2
        %contstruct.firstpoint = get(gca,'currentpoint');
        %contstruct.startx=contstruct.firstpoint(1,1);
		%contstruct.starty=contstruct.firstpoint(1,2);
		%set(findobj(gca,'tag','contdata'),'userdata',contstruct);
    elseif contstruct.startpoint == 2
        contstruct.firstpoint = get(gca,'currentpoint');
        contstruct.startpoint= 1;
        contstruct.startx=contstruct.firstpoint(1,1);
		contstruct.starty=contstruct.firstpoint(1,2);
		set(findobj(gca,'tag','contdata'),'userdata',contstruct);
    else 
   		contstruct.secondpoint = get(gca,'currentpoint');
		contstruct.endx=contstruct.secondpoint(1,1);
		contstruct.endy=contstruct.secondpoint(1,2);
        
        
%		fprintf('X:   %0.2f   Y:   %0.2f\n',-(startx-endx),starty-endy);
		contstruct.wx=contstruct.windowval+(contstruct.starty-contstruct.endy).*contstruct.scaley;
		contstruct.cx=contstruct.centerval+(contstruct.startx-contstruct.endx).*contstruct.scalex;
        
		
		if contstruct.cx < 0
			contstruct.cx = 0.0001;
		end

		if contstruct.wx < 0
			contstruct.wx = 0.0001;
		end
		bot=round(contstruct.cx-contstruct.wx/2);
		top=round(contstruct.cx+contstruct.wx/2);

		if bot < 0
			bot =0;
		end		
        
        clim=contstruct.clim(1)+(contstruct.clim(2)-contstruct.clim(1))*[contstruct.cx-contstruct.wx/2 contstruct.cx+contstruct.wx/2];
        %set(gca,'clim',[contstruct.cx-contstruct.wx/2 contstruct.cx+contstruct.wx/2]);		
		%col=findobj(gcf,'type','axes')
        %d=get(gcf,'children');
        %[b]=size(d)
        %clim=get(gca,'clim')
        %lim2=get(a(2),'clim')
        if clim(1) == clim(2), clim=[clim(1) clim(2)+0.0001];end
        caxis([clim]);
   
        hbar=cmap(gca);
        %[sz y]=size(get(gcf,'colormap'));
        %b=(sz-1)./(contstruct.clim(2)-contstruct.clim(1));
        %a=1-b*contstruct.clim(1)
        %climbar=caxis(gca)
        %climbar=a+b.*clim
        %set(d(2),'ylim',climbar)
        
        if clim(2) > contstruct.clim(2),clim(2)=contstruct.clim(2);end
        if clim(1) < contstruct.clim(1),clim(1)=contstruct.clim(1);end
        if clim(1) > clim(2),return,end
        set(hbar,'ylim',clim)
        
        set(findobj(gca,'tag','contdata'),'userdata',contstruct);
		%set(findobj(gca,'tag','contdata'),'string',sprintf('C: %d W: %d',round(contstruct.cx),round(contstruct.wx)));
	    end
      
    case 'stop'
   	set(gcf,'WindowButtonMotionFcn','')
   	set(gcf,'WindowButtonUpFcn','')
    set(gcbf,'pointer','arrow')
	contstruct=get(findobj(gca,'tag','contdata'),'userdata');
   	contstruct.startpoint=1;
   	contstruct.centerval=contstruct.cx;
   	contstruct.windowval=contstruct.wx;
	set(findobj(gca,'tag','contdata'),'userdata',contstruct);
	set(findobj(gca,'tag','contdata'),'string','');
    
    clim=contstruct.clim(1)+(contstruct.clim(2)-contstruct.clim(1))*[contstruct.cx-contstruct.wx/2 contstruct.cx+contstruct.wx/2];
		
    if clim(2) > contstruct.clim(2),clim(2)=contstruct.clim(2);end
    if clim(1) < contstruct.clim(1),clim(1)=contstruct.clim(1);end
  
   
        
    %hbar=colorbar;set(hbar,'ylim',clim);
    %lim=get(hbar,'clim')
    
    case 'reset' 
        
        contstruct=get(findobj(gca,'tag','contdata'),'userdata');
       
        set(gca,'clim',contstruct.clim);
        caxis([contstruct.clim]);
        contstruct.windowval=0.5;
        contstruct.centerval=0.25;
        contstruct.cx=0.5;
        contstruct.wx=0.5;
        contstruct.startpoint=2;
        clear contstruct.endx contstruct.endy contstruct.startx contstruct.starty 
        clear contstruct.firstpoint contstruct.secondpoint contstruct.lastpoint
        set(findobj(gca,'tag','contdata'),'userdata',contstruct);
      
        hbar=cmap(gca);
       
        %contstruct.startpoint == 1
        set(hbar,'clim',contstruct.clim);
        set(hbar,'ylim',contstruct.clim);    
        
end
   