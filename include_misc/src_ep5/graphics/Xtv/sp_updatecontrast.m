function sp_updatecontrast(action)
%global windowval centerval startpoint startx starty imsize currentgamma cx wx
switch action

case 'start'
    mouse=get(gcf,'selectiontype');
    if strcmp('normal',mouse)    
    set(gcbf,'WindowButtonMotionFcn','sp_updatecontrast move')
   	set(gcbf,'WindowButtonUpFcn','sp_updatecontrast stop')
	set(gcbf,'pointer','fleur')
    end
    if strcmp('alt',mouse)
    sp_updatecontrast reset
    end    
    if strcmp('extend',mouse)
    set(gcbf,'WindowButtonUpFcn','')
    %hbar=colorbar;set(hbar,'ylim',contstruct.clim);    
    end
    
    sp_updatecontrast move

case 'move'
	contstruct=get(findobj(gca,'tag','contdata'),'userdata');
    
        if isempty( contstruct )
            return
        end
         
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
        
        if isequal(clim(1),clim(2)) ==0
        caxis([clim]);
        end   
        handles=guidata(gcf);
        hbar=sp_cbar(handles);
        %hbar=cmap(gca);
        %[sz y]=size(get(gcf,'colormap'));
        %b=(sz-1)./(contstruct.clim(2)-contstruct.clim(1));
        %a=1-b*contstruct.clim(1)
        %climbar=caxis(gca)
        %climbar=a+b.*clim
        %set(d(2),'ylim',climbar)
        
        if clim(2) > contstruct.clim(2),clim(2)=contstruct.clim(2);end
        if clim(1) < contstruct.clim(1),clim(1)=contstruct.clim(1);end
        if clim(1) > clim(2),return,end
       
        if isequal(clim(1),clim(2)),
        set(hbar,'ylimmode','auto');    
        else
        set(hbar,'ylim',clim)
        end
        
        set(findobj(gca,'tag','contdata'),'userdata',contstruct);
		%set(findobj(gca,'tag','contdata'),'string',sprintf('C: %d W: %d',round(contstruct.cx),round(contstruct.wx)));
	    end
      
    case 'stop'
   	set(gcf,'WindowButtonMotionFcn','')
   	set(gcf,'WindowButtonUpFcn','')
    set(gcbf,'pointer','arrow')
	contstruct=get(findobj(gca,'tag','contdata'),'userdata');
   	contstruct.startpoint=1;
        if isfield( contstruct, 'cx' )
   	    contstruct.centerval=contstruct.cx;
        end
        if isfield( contstruct, 'wx' )
            contstruct.windowval=contstruct.wx;
        end
	set(findobj(gca,'tag','contdata'),'userdata',contstruct);
	set(findobj(gca,'tag','contdata'),'string','');
    
   
    if isfield( contstruct, 'clim' )
        clim=contstruct.clim(1)+(contstruct.clim(2)-contstruct.clim(1))*[contstruct.cx-contstruct.wx/2 contstruct.cx+contstruct.wx/2];
    else
        return
    end
		
    if clim(2) > contstruct.clim(2),clim(2)=contstruct.clim(2);end
    if clim(1) < contstruct.clim(1),clim(1)=contstruct.clim(1);end
  
    %%%%%%%%%%%%%%%%%%%%%%%%
    % for scanp
    
    impix(gcf) 
    handles=guidata(gcf);
    
    
        
     minval=clim(1);
     maxval=clim(2);
      
     set(handles.min_edit,'String',strrep(sprintf('%7.3g',minval),'',''));
     set(handles.max_edit,'String',strrep(sprintf('%7.3g',maxval),'',''));
%     
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %hbar=colorbar;set(hbar,'ylim',clim);
    %lim=get(hbar,'clim')
    
    case 'reset' 
        
        contstruct=get(findobj(gca,'tag','contdata'),'userdata');
        
        if isequal(contstruct.clim(1),contstruct.clim(2)),
        else
        set(gca,'clim',contstruct.clim);
        
        
        handles=guidata(gcf);
        caxis([handles.minval handles.maxval]);
        
        
        contstruct.windowval=1;
        contstruct.centerval=0.5;
        contstruct.cx=0.5;
        contstruct.wx=1;
        contstruct.startpoint=2;
        clear contstruct.endx contstruct.endy contstruct.startx contstruct.starty 
        clear contstruct.firstpoint contstruct.secondpoint contstruct.lastpoint
        set(findobj(gca,'tag','contdata'),'userdata',contstruct);
        handles=guidata(gcf);
        hbar=sp_cbar(handles);
        set(hbar,'clim',contstruct.clim);
        set(hbar,'ylim',contstruct.clim);  
             
        set(handles.min_edit,'String',strrep(sprintf('%7.3g',contstruct.clim(1)),'',''));            
        set(handles.max_edit,'String',strrep(sprintf('%7.3g',contstruct.clim(2)),'',''));
        set(findobj(gca,'tag','contdata'),'userdata',contstruct);
    
        %contstruct.startpoint == 1
        set(hbar,'clim',contstruct.clim);
        set(hbar,'ylim',contstruct.clim);    
        end
end
   
