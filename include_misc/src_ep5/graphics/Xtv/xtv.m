function varargout= xtv(varargin)
% xtv Application M-file for xtv.fig
%    FIG = xtv launch xtv GUI.
%    xtv('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 28-Jan-2013 15:51:40

if nargin==0, varargin{1}=[];data=[];end
if isnumeric(varargin{1})  % LAUNCH GUI
    %disp('bitte eine 2D oder 3D Matrix der Form (NY x NX x NZ)
    %übergeben!')
    data = squeeze(varargin{1});
if nargin==2, dim=varargin{2};else dim=1;end 
    fig = openfig(mfilename,'new');
    	%set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
    
	% Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);
    %handles.orgdata = data./max(abs(data(:)));
    handles.orgdata = data;
    handles.datasize = size(handles.orgdata);
    handles.dim     = dim;
    handles.flip    = [0 0 0];
    handles.perm    = [1 2 3];
    handles.FOV     = [200 200 200];
    load pointer
    handles.pointer=pointer;
    %handles.z = ;
    handles.xi=[];
    handles.yi=[];
    
    handles.maxval = max(abs(handles.orgdata(:)));
    handles.minval = min(abs(handles.orgdata(:)));
    
    guidata(fig, handles);
    %ishandle(handles.slider1)
    %set(handles.slider1,'visible','off')
    %set(handles.figure1,'CurrentAxes',handles.axes1)
    %guidata(fig,handles);
    %impix(handles.figure1)
    %handles=getdata(handles);
    %guidata(fig,handles);
    %if isempty(lasterr)==1, return, end
    pushbutton10_Callback(handles.pushbutton10,[],handles)
 
%------------------------------------------------------------------------------    
%   So können die Menü-Funktionen beispielsweise durch CTRL+C aufgerufen werden
    set(handles.print_submenu,'Accelerator','p')
    set(handles.close_submenu,'Accelerator','c')
    set(handles.write_submenu,'Accelerator','w')
    set(handles.load_submenu,'Accelerator','l')

%------------------------------------------------------------------------------
%   Diese Funktion ist das Herz des ganzen 
   %set(handles.figure1,'HandleVisibility','on')
   %impix(handles.figure1)

%------------------------------------------------------------------------------    
    
    if nargout > 0
		varargout{1} = fig;
	end
end

if ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end 
end



%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switcset(handles.figure1,'CurrentAxes',handles.axes1hyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.



% --------------------------------------------------------------------
function varargout = imgnr_slider_Callback(h, eventdata, handles, varargin)



% We need no reset
set(handles.ToggleReset,'Value',0);


val=get(handles.imgnr_slider,'Value');
imgnr=round(1+(handles.datasize(3)-1)*val);
imgnrstr=strrep(sprintf('%d',imgnr),' ','');
set(handles.imgnr_edit,'String',imgnrstr)
newval=(imgnr-1)/(handles.datasize(3)-1);
set(handles.imgnr_slider,'Value',newval)
handles.z=imgnr;
set(handles.figure1,'CurrentAxes',handles.axes1);
%img=handles.data(:,:,round(imgnr));
guidata(h,handles)
%minval=str2double(get(handles.min_edit,'String'));
%maxval=str2double(get(handles.max_edit,'String'));
%img(find(img<minval))=minval;
%img(find(img>maxval))=maxval;
%handles.img=img;
set(handles.figure1,'HandleVisibility','on')
%imagesc(img)

update_image(handles,0)
%if ishandle(get(handles.axes5,'Children')), colorbar(handles.axes5),end

if isempty(get(handles.figure1,'WindowButtonMotionFcn')),
    if isempty(handles.xi)==0 & isempty(handles.yi)==0
    [handles.roi]=roipoly(handles.datasize(1),handles.datasize(2),handles.xi,handles.yi);
    line(handles.xi,handles.yi,'Color','w')
    plotroi(handles)
        elseif isempty(handles.x1)==0 & isempty(handles.y1)==0 &...
            isempty(handles.x2)==0 & isempty(handles.y2)==0 &...
            isempty(handles.x3)==0 & isempty(handles.y3)==0 &...
            isempty(handles.x4)==0 & isempty(handles.y4)==0 
        calcsnr(handles)
    end
end   
%update_image(handles,0)
%set(handles.axes1,'UserData',handles.data);
%set(handles.figure1,'HandleVisibility','off')
guidata(h,handles)

%----------------------------------------------------------------------

function varargout = imgnr_edit_Callback(h, eventdata, handles, varargin)

% We need no reset
set(handles.ToggleReset,'Value',0);


imgnr=str2double(get(h,'String'));
handles.z=imgnr;
val=(imgnr-1)/(handles.datasize(3)-1);
set(handles.imgnr_slider,'Value',val)
set(handles.figure1,'CurrentAxes',handles.axes1);

update_image(handles,0)
set(handles.axes1,'UserData',handles.data);
guidata(h,handles)

% --------------------------------------------------------------------
 
% Radio buttons

% --------------------------------------------------------------------

function varargout = real_radio_Callback(h, eventdata, handles, varargin)
off=[handles.abs_radio,handles.imag_radio,handles.phase_radio];
mutual_exclude(off)
set(handles.real_radio,'enable','off')

% we need a reset
set(handles.ToggleReset,'Value',1);


handles.data = real(permute(handles.orgdata,handles.perm));

if handles.flip(1)
    handles.data = flipdim(handles.data,1);
end

if handles.flip(2)
    handles.data = flipdim(handles.data,2);
end

if handles.flip(3)
    handles.data = flipdim(handles.data,3);
end

handles.maxval=max(handles.data(:));
handles.minval=min(handles.data(:));
minvalstr=strrep(sprintf('%7.3g',handles.minval),'','');
maxvalstr=strrep(sprintf('%7.3g',handles.maxval),'','');

caxis([handles.minval handles.maxval]);

set(handles.min_edit,'String',minvalstr)
set(handles.max_edit,'String',maxvalstr)


state=get(handles.dreiD_toggle,'Value');
update_image(handles,state)
guidata(h,handles)  



% --------------------------------------------------------------------
function varargout = abs_radio_Callback(h, eventdata, handles, varargin)

off=[handles.real_radio,handles.imag_radio,handles.phase_radio];
mutual_exclude(off)
set(handles.abs_radio,'enable','off')
set(handles.abs_radio,'value',1)

% we need a reset
set(handles.ToggleReset,'Value',1);

handles.data = abs(permute(handles.orgdata,handles.perm));

if handles.flip(1)
    handles.data = flipdim(handles.data,1);
end

if handles.flip(2)
    handles.data = flipdim(handles.data,2);
end

if handles.flip(3)
    handles.data = flipdim(handles.data,3);
end

handles.maxval=max(handles.data(:));
handles.minval=min(handles.data(:));

caxis([handles.minval handles.maxval]);

minvalstr=strrep(sprintf('%7.3g',handles.minval),'','');
maxvalstr=strrep(sprintf('%7.3g',handles.maxval),'','');
set(handles.min_edit,'String',minvalstr)
set(handles.max_edit,'String',maxvalstr)

state=get(handles.dreiD_toggle,'Value');
update_image(handles,state)
guidata(h,handles)  

%handles.state='abs';

% if ndims(handles.orgdata) == 3
% perm = handles.perm;
% handles.data = abs(permute(handles.orgdata,perm));
% else
%     handles.data = handles.orgdata
% end
%set(handles.figure1,'CurrentAxes',handles.axes1)
%caxis= [0 1];

%data = abs(handles.data);
% minval=min(data(:)); minvalstr=strrep(sprintf('%7.3g',minval),'','');
% maxval=max(data(:)); maxvalstr=strrep(sprintf('%7.3g',maxval),'','');
% handles.minval=minval;
% handles.maxval=maxval;
% set(handles.min_edit,'String',minvalstr)
% set(handles.max_edit,'String',maxvalstr)
%handles.data = abs(squeeze(handles.orgdata));
%handles.maxval=max(handles.data(:));
%handles.minval=min(handles.data(:));
%state=get(handles.dreiD_toggle,'Value');
%set(handles.axes1,'UserData',abs(handles.data));
%update_image(handles,state)
%guidata(h,handles)    


% --------------------------------------------------------------------

function varargout = imag_radio_Callback(h, eventdata, handles, varargin)
off=[handles.real_radio,handles.abs_radio,handles.phase_radio];
mutual_exclude(off)
set(handles.imag_radio,'enable','off')

% we need a reset
set(handles.ToggleReset,'Value',1);


handles.data = imag(permute(handles.orgdata,handles.perm));

if handles.flip(1)
    handles.data = flipdim(handles.data,1);
end

if handles.flip(2)
    handles.data = flipdim(handles.data,2);
end

if handles.flip(3)
    handles.data = flipdim(handles.data,3);
end

handles.maxval=max(handles.data(:));
handles.minval=min(handles.data(:));

caxis([handles.minval handles.maxval]);

minvalstr=strrep(sprintf('%7.3g',handles.minval),'','');
maxvalstr=strrep(sprintf('%7.3g',handles.maxval),'','');
set(handles.min_edit,'String',minvalstr)
set(handles.max_edit,'String',maxvalstr)

state=get(handles.dreiD_toggle,'Value');
update_image(handles,state)
guidata(h,handles)  



% --------------------------------------------------------------------
function varargout = phase_radio_Callback(h, eventdata, handles, varargin)
off=[handles.real_radio,handles.abs_radio,handles.imag_radio];
mutual_exclude(off)
set(handles.phase_radio,'enable','off')

% we need a reset
set(handles.ToggleReset,'Value',1);


handles.data = angle(permute(handles.orgdata,handles.perm));

if handles.flip(1)
    handles.data = flipdim(handles.data,1);
end

if handles.flip(2)
    handles.data = flipdim(handles.data,2);
end

if handles.flip(3)
    handles.data = flipdim(handles.data,3);
end

handles.maxval = max(handles.data(:));
handles.minval = min(handles.data(:));

caxis([handles.minval handles.maxval]);


minvalstr=strrep(sprintf('%7.3g',handles.minval),'','');
maxvalstr=strrep(sprintf('%7.3g',handles.maxval),'','');
set(handles.min_edit,'String',minvalstr)
set(handles.max_edit,'String',maxvalstr)

state=get(handles.dreiD_toggle,'Value');
update_image(handles,state)
guidata(h,handles)  



%---------------------------------------------------------------------
function mutual_exclude(off)
set(off,'Value',0)
set(off,'enable','on')
% --------------------------------------------------------------------


% --------------------------------------------------------------------


function varargout = dreiD_toggle_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.dreiD_toggle.
set(handles.figure1,'CurrentAxes',handles.axes1)
state=get(h,'Value');
if state==1, 
set(handles.imgnr_slider,'Visible','off')
set(handles.imgnr_edit,'Visible','off')
set(handles.roi_submenu,'enable','off')
set(handles.snr_submenu,'enable','off')
else
set(handles.imgnr_slider,'Visible','on')
set(handles.imgnr_edit,'Visible','on')
set(handles.roi_submenu,'enable','on')
set(handles.snr_submenu,'enable','on')
end

set(handles.ToggleReset,'Value',1);
set(handles.ZoomToggle,'Value',0);
guidata(h,handles)
update_image(handles,state)
guidata(h,handles)


% --------------------------------------------------------------------


% --------------------------------------------------------------------
% 
%   File-Menu
%
% --------------------------------------------------------------------
function varargout = file_menu_Callback(h, eventdata, handles, varargin)
if isempty(get(handles.figure1,'WindowButtonMotionFcn'))
    set(handles.close_submenu,'Enable','off')
else
    set(handles.close_submenu,'Enable','on')
end

if isempty(get(handles.axes1,'Children'))
    set(handles.print_submenu,'Enable','off')
else
    set(handles.print_submenu,'Enable','on')
end

% --------------------------------------------------------------------

function varargout = print_submenu_Callback(h, eventdata, handles, varargin)
print -f handles.figure1 -noui -v 

% --------------------------------------------------------------------

function varargout = close_submenu_Callback(h, eventdata, handles, varargin)
delete(handles.figure1)

% --------------------------------------------------------------------
function varargout = write_submenu_Callback(h, eventdata, handles, varargin)
set(handles.figure1, 'WindowButtonMotionFcn', '', ...
    'Interruptible', 'on', 'BusyAction', 'Queue');
[fname,pathname]=uiputfile('*.tif;*.ps;*.eps;*.jpg');
        a=findstr(fname,'.');
 if isempty(a), format='tif',else format=fname(a(end)+1:end);end
 full_path=[pathname fname];
 if strcmp(format,'tif') | strcmp(format,'tiff') 
        print(handles.figure1,'-dtiff',full_path)
 elseif  strcmp(format,'ps')
        print(handles.figure1,'-dpsc2',full_path)
 elseif  strcmp(format,'eps')
        print(handles.figure1,'-depsc2',full_path)
 elseif  strcmp(format,'jpg')
        print(handles.figure1,'-djpeg',full_path)
 else 
     uiwait(msgbox('Kein gültiges Grafikformat!','Write','error','modal'))
 end
    impix(handles.figure1)
% --------------------------------------------------------------------
function varargout = load_submenu_Callback(h, eventdata, handles, varargin)
    [fname,pathname]=uigetfile('*.mat');
        full_path=[pathname fname];
        data=fload(full_path);
        handles.orgdata=data;
        handles.data=abs(data);
        guidata(h,handles)
        handles=getdata(handles);
        impix(handles.figure1)
        guidata(h,handles)

        
% --------------------------------------------------------------------
%
%        Tools-Menu
%
% --------------------------------------------------------------------

function varargout = tools_menu_Callback(h, eventdata, handles, varargin)

if get(handles.abs_radio,'Value')==0 
    set(handles.snr_submenu,'Enable','off')
else
    set(handles.snr_submenu,'Enable','on')
end

% --------------------------------------------------------------------
function varargout = colorbar_submenu_Callback(h, eventdata, handles, varargin)
a=get(handles.axes5,'Children');
state=get(handles.axes5,'Visible');
maxvalstr=get(handles.max_edit,'String');maxval=str2double(maxvalstr);
minvalstr=get(handles.min_edit,'String');minval=str2double(minvalstr);
if isequal('off',state), set(handles.axes5,'Visible','on'),...
        handles.axes5=sp_cbar(handles);
        %colorbar(handles.axes5,'peer',handles.axes1),...
        %set(handles.axes5,'ylim',[minval maxval]),
end
if isequal('on',state), set(handles.axes5,'Visible','off'),delete(a),end

% --------------------------------------------------------------------
function varargout = pixval_submenu_Callback(h, eventdata, handles, varargin)
pixval(handles.figure1)
% --------------------------------------------------------------------
function varargout = xadjust_submenu_Callback(h, eventdata, handles, varargin)
xadjust
% --------------------------------------------------------------------

function varargout = roi_submenu_Callback(h, eventdata, handles, varargin)
set(handles.figure1, 'WindowButtonMotionFcn', '', ...
    'Interruptible', 'off', 'BusyAction', 'Queue');
set(handles.figure1,'CurrentAxes',handles.axes1)
[handles.roi,xi,yi]=roipoly;
line(xi,yi,'Color','w')
handles.xi=xi;
handles.yi=yi;
guidata(h,handles)
plotroi(handles)
set(handles.figure1,'Pointer','arrow')
wait(handles)
%handles.xi=[];
%handles.yi=[];
guidata(h,handles)



function varargout = snr_submenu_Callback(h, eventdata, handles, varargin)

if isfield(handles,'x1')==0 & isfield(handles,'y1')==0 ...
        isfield(handles,'x2')==0 & isfield(handles,'y2')==0 & ...
        isfield(handles,'x3')==0 & isfield(handles,'y3')==0
    uiwait(msgbox('Bitte 3 ROIs im Rauschen festlegen','SNR','modal'))
    set(handles.figure1,'CurrentAxes',handles.axes1)
    [roi1,x1,y1]=roipoly; noise1=handles.data(:,:,handles.z).*double(roi1); 
        noise(1)=std(noise1(find(noise1~=0)));
        handles.x1=x1;
        handles.y1=y1;
    [roi2,x2,y2]=roipoly; noise2=handles.data(:,:,handles.z).*double(roi2); 
        noise(2)=std(noise2(find(noise2~=0)));
        handles.x2=x2;
        handles.y2=y2;
    [roi3,x3,y3]=roipoly; noise3=handles.data(:,:,handles.z).*double(roi3); 
        noise(3)=std(noise3(find(noise3~=0)));
    handles.x3=x3;
    handles.y3=y3;
   
else 
[roi1]=double(roipoly(handles.datasize(1),handles.datasize(2),handles.x1,handles.y1)); 
        noise1=handles.data(:,:,handles.z).*roi1; noise(1)=std(noise1(find(noise1~=0)));
[roi2]=double(roipoly(handles.datasize(1),handles.datasize(2),handles.x2,handles.y2)); 
        noise2=handles.data(:,:,handles.z).*roi2; noise(2)=std(noise2(find(noise2~=0)));
[roi3]=double(roipoly(handles.datasize(1),handles.datasize(2),handles.x3,handles.y3)); 
        noise3=handles.data(:,:,handles.z).*roi3; noise(3)=std(noise3(find(noise3~=0)));
end
noise=mean(noise(:));
guidata(h,handles)
uiwait(msgbox('Bitte 1 ROI im Signal festlegen','SNR','modal'))
set(handles.figure1,'CurrentAxes',handles.axes1)
[roi4,x4,y4]=roipoly; signal=handles.data(:,:,handles.z).*double(roi4); 
    line(x4,y4,'Color','y')
    signal=mean(signal((find(signal~=0))));
    handles.x4=x4;
    handles.y4=y4;
    handles.xi=[];
    handles.yi=[];
guidata(h,handles)
SNR=0.655*(signal/noise);
str = sprintf('Nr.= %g, SNR= %6.2f ',handles.z,SNR);
set(handles.update,'String',str)
set(handles.figure1, 'WindowButtonMotionFcn', '', ...
    'Interruptible', 'on', 'BusyAction', 'Queue');
set(handles.figure1,'Pointer','arrow')
wait(handles)


function varargout = zoom_submenu_Callback(h, eventdata, handles, varargin)

zoomout(handles)
handles.xLim = get(handles.axes1,'Xlim');
handles.yLim = get(handles.axes1,'ylim');
handles.zoom = 1;
guidata(h,handles)

function varargout = distance_submenu_Callback(h, eventdata, handles, varargin)

kids=get(handles.axes1,'children');
b=get(kids(:),'type');
a=strcmp(b,'line');
for k=1:size(b)-1, delete(kids(k));end
dist=calcdist(handles.axes1);
sprintf('Distance= %g',dist)
set(handles.update,'String',sprintf('Distance= %g',dist))
set(handles.figure1,'WindowButtonMotionFcn','')
pause(2)

update_image(handles,0)


function varargout = movie_submenu_Callback(h, eventdata, handles, varargin)

set(handles.figure1,'CurrentAxes',handles.axes1)

currentframe = handles.z;
cmap=colormap;
%mov=immovie(length(cmap).*mat2gray(handles.data),cmap);
%update_image(handles,0)
prompt={'Enter frame rate:',...
        'Enter the filename:'};
name='Movie Options';
numlines=1;
defaultanswer={'20','movie.avi'};
answer = inputdlg(prompt,name,numlines,defaultanswer);

%FrameRate = inputdlg('Choose FrameRate. Default is 20','FrameRate',2)

vidObj = VideoWriter(answer{2});
vidObj.FrameRate = str2double(answer{1});

open(vidObj);

for f = 1:handles.datasize(3)
    handles.z = f;
    update_image(handles,0)
    currFrame = getframe;
    writeVideo(vidObj,currFrame);
end

 close(vidObj);

handles.z = currentframe;
update_image(handles,0)



% --------------------------------------------------------------------
function MIP_submenu_Callback(h, eventdata, handles, varargin)
% hObject    handle to MIP_submenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt={'Enter first Frame:',...
        'Enter last Frame:'};
name='MIP Options';
numlines=1;
defaultanswer={'1',num2str(handles.datasize(3))};
answer = inputdlg(prompt,name,numlines,defaultanswer);
try
    ima = max(handles.data(:,:,str2double(answer{1}):str2double(answer{2})),[],3);
catch ME
    error('Input was not correct, choose appropriate values');
end
imghandle=tv(ima,3,[],1);
%guidata(h,handles)



% --------------------------------------------------------------------
% 
%   Help-Menu
%
% --------------------------------------------------------------------

function varargout = help_menu_Callback(h, eventdata, handles, varargin)

function varargout = info_submenu_Callback(h, eventdata, handles, varargin)
message=[char(169) '  by Felix Breuer 2002']
uiwait(msgbox(message,'Info','modal'))

function varargout = snrhelp_submenu_Callback(h, eventdata, handles, varargin)
CreateStruct.WindowStyle='replace';
CreateStruct.Interpreter='tex';
uiwait(msgbox('{\itSNR}=0.655 \cdot {<\itS>}/{<\sigma_R>}','SNR-Info','help',CreateStruct));

% --------------------------------------------------------------------
% --------------------------------------------------------------------
function plotroi(handles)
z=handles.z;


if size(handles.datasize,2)==2, total=1;else total=handles.datasize(3);end




for k=1:total, handles.roidata(:,:,k)=double(handles.roi).*squeeze(handles.data(:,:,k));end


set(handles.figure1,'CurrentAxes',handles.axes2)
set(handles.axes2,'XlimMode','auto')
set(handles.axes2,'YlimMode','auto')

currdata=handles.roidata(:,:,z);
hist(currdata(currdata>0),100),set(handles.text9,'string','histogram roi')
set(handles.figure1,'CurrentAxes',handles.axes3)
set(handles.axes2,'XlimMode','auto')
set(handles.axes2,'XtickMode','auto')




prodata=reshape(handles.roidata,handles.datasize(1)*handles.datasize(2),total);
for k=1:total, tmp=prodata(:,k); tmp=tmp(tmp~=0); meanpro(k)=mean(tmp);end
for k=1:total, tmp=prodata(:,k); tmp=tmp(tmp~=0); stdpro(k)=std(tmp);end

if total>1, 
    plot(meanpro,'kx:'), set(handles.text10,'string','z profile (mean)')
    stdval=std(squeeze(meanpro(:)));
    maxval=max(squeeze(meanpro(:)));
    minval=min(squeeze(meanpro(:)));
    set(handles.axes3,'Ylim',[minval-stdval maxval+stdval]);

    step=round(handles.datasize(3)/10);
    set(handles.axes3,'Xlim',[1 handles.datasize(3)],'Xtick',[1:step:handles.datasize(3)])

    ylim=get(handles.axes3,'Ylim');
    line([z z],[ylim(1) ylim(2)])

    set(handles.figure1,'CurrentAxes',handles.axes4)

    plot(stdpro,'kx:'),set(handles.text11,'string','z profile (std)')
    
    stdval=std(squeeze(stdpro(:)));
    maxval=max(squeeze(stdpro(:)));
    minval=min(squeeze(stdpro(:)));
    set(handles.axes4,'Ylim',[minval-stdval maxval+stdval]);
    step=round(handles.datasize(3)/10);
    set(handles.axes4,'Xlim',[1 handles.datasize(3)],'Xtick',[1:step:handles.datasize(3)])

    ylim=get(handles.axes4,'Ylim');
    line([z z],[ylim(1) ylim(2)])
end

str = sprintf('z= %g, points= %g,  %6.4f %s %6.4f',...
    z,size(handles.roi(handles.roi~=0),1),meanpro(z),char(49+128),stdpro(z));
set(handles.update,'String',str)

set(handles.figure1, 'WindowButtonMotionFcn', '', ...
    'Interruptible', 'on', 'BusyAction', 'Queue');
%update_image(handles,0)



function wait(handles)
key=get(handles.figure1,'CurrentCharacter');

waitfor(handles.imgnr_slider,'value')

%a
%waitfor(handles.fgure1,'SelectionType')
%w=waitforbuttonpress
% if w==1
%    
%     if strcmp(key,'p'), print -f handles.figure1 -v, w=waitforbuttonpress, end
%     if strcmp(key,'w')
%         [fname,pathname]=uiputfile('*.tif','*.tiff');
%         full_path=[pathname fname];
%         print(handles.figure1,'-dtiff',full_path)
%         w=waitforbuttonpress;
%     end    
% end
set(handles.figure1,'CurrentAxes',handles.axes1)
set(handles.figure1,'SelectionType','normal')
update_image(handles,0)

%impix(handles.figure1)
%set(handles.figure1,'HandleVisibility','off')



function handles = getdata(handles)

handles.data = abs(squeeze(handles.orgdata));

if ndims(handles.data)==3;
    perm=circshift([1 2 3],[0 -handles.dim]);
    handles.data = permute(handles.data,perm);
    handles.datasize = size(handles.data);
    %handles.FOV      = size(handles.data);
    handles.perm     = perm;
    set(handles.pushbutton_Permute,'visible','on')
    set(handles.pushbutton_flipBF,'visible','on')
    set(handles.imgnr_slider,'Visible','on')
    set(handles.imgnr_edit,'visible','on')
    set(handles.dreiD_toggle,'Visible','on')
    set(handles.frame9,'Visible','on')
    set(handles.imgnr_slider,'SliderStep', ...
    [1./(handles.datasize(3)-1) 1./(handles.datasize(3)-1)])
    
    handles.z = floor(handles.datasize(3)/2)+1;
    val=(handles.z-1)/(handles.datasize(3)-1);

    set(handles.imgnr_slider,'Value',val)
    imgnrstr=strrep(sprintf('%d',handles.z),' ','');
    set(handles.imgnr_edit,'String',imgnrstr)
            
end

if ndims(handles.data)==2;
    handles.datasize = [size(squeeze(handles.data)) 1];
    handles.FOV(3)   = 1;
    handles.z = 1;
    set(handles.pushbutton_Permute,'visible','off')
    set(handles.pushbutton_flipBF,'visible','off')
    
    set(handles.imgnr_edit,'visible','off')
    set(handles.imgnr_slider,'Visible','off')
    set(handles.dreiD_toggle,'Visible','off')
    set(handles.frame9,'Visible','off')
end

set(handles.ToggleReset,'visible','off')

map=get(0, 'DefaultFigurecolormap');
colormap(map);

minval=min(abs(handles.data(:))); minvalstr=strrep(sprintf('%7.3g',minval),'','');
maxval=max(abs(handles.data(:))); maxvalstr=strrep(sprintf('%7.3g',maxval),'','');
% 
handles.minval=minval;
handles.maxval=maxval;
%update_image(handles,0)

set(handles.min_edit,'String',minvalstr)
set(handles.max_edit,'String',maxvalstr)

ima=handles.data(:,:,handles.z);
imghandle=tv(ima,3,[],1);

handles.xLim = get(handles.axes1,'Xlim');
handles.yLim = get(handles.axes1,'Ylim');

    
abs_radio_Callback(handles.abs_radio,[],handles)



%set(handles.axes1,'UserData',handles.data);

%set(handles.axes1,'HandleVisibility','Callback')

function calcsnr(handles)

set(handles.figure1,'CurrentAxes',handles.axes1)

[roi1]=double(roipoly(handles.datasize(1),handles.datasize(2),handles.x1,handles.y1)); 

noise1=handles.data(:,:,handles.z).*roi1; noise(1)=std(noise1(find(noise1~=0)));

[roi2]=double(roipoly(handles.datasize(1),handles.datasize(2),handles.x2,handles.y2)); 
        noise2=handles.data(:,:,handles.z).*roi2; noise(2)=std(noise2(find(noise2~=0)));
[roi3]=double(roipoly(handles.datasize(1),handles.datasize(2),handles.x3,handles.y3)); 
        noise3=handles.data(:,:,handles.z).*roi3; noise(3)=std(noise3(find(noise3~=0)));
noise=mean(noise(:));
%msgbox('Bitte 1 ROI im Signal festlegen')
[roi4]=double(roipoly(handles.datasize(1),handles.datasize(2),handles.x4,handles.y4)); 
signal=handles.data(:,:,handles.z).*roi4;size(signal); signal=mean(signal(find(signal~=0)));
line(handles.x4,handles.y4,'Color','y')
SNR=0.655*(signal/noise);
str = sprintf('Nr.= %g, SNR= %6.2f ',handles.z,SNR);
set(handles.update,'String',str)


% --------------------------------------------------------------------
function varargout = pushbutton10_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton10.
%disp('pushbutton10 Callback not implemented yet.')
set(handles.figure1,'CurrentAxes',handles.axes1)
handles=getdata(handles);
%caxis([0 1])
caxis([handles.minval handles.maxval]);
impix(handles.figure1)
guidata(h,handles);


% --------------------------------------------------------------------
function varargout = min_edit_Callback(h, eventdata, handles, varargin)

maxvalstr=get(handles.max_edit,'String');maxval=str2double(maxvalstr);
minvalstr=get(handles.min_edit,'String');minval=str2double(minvalstr);

%minimum=handles.minval;
%maximum=handles.maxval;

if isnan(maxval) | isnan(minval) | ... 
   minval>maxval 
    %errordlg('You must enter a numeric value','Bad Input','modal')
    set(handles.min_edit,'String',strrep(sprintf('%7.3g',handles.oldmin),'',''));
    return
end

%if minval<minimum, minval=minimum; ...
    set(handles.min_edit,'String',strrep(sprintf('%7.3g',minval),'',''));
%end
%if maxval>maximum, maxval=maximum; ... 
    set(handles.max_edit,'String',strrep(sprintf('%7.3g',maxval),'',''));
%end

caxis([minval maxval]./(max(handles.data(:))-min(handles.data(:))));

handles.oldmax=maxval;
handles.oldmin=minval;

set(handles.figure1,'CurrentAxes',handles.axes1);
state=get(handles.dreiD_toggle,'value');
guidata(h,handles)
%scale_image(handles,state)
update_image(handles,state)
handles.axes5=sp_cbar(handles);
guidata(h,handles)

% --------------------------------------------------------------------
function varargout = max_edit_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.max_edit.
%disp('max_edit Callback not implemented yet.')
minvalstr=get(handles.min_edit,'String');minval=str2double(minvalstr);
maxvalstr=get(handles.max_edit,'String');maxval=str2double(maxvalstr);
%minimum=handles.minval;
%maximum=handles.maxval;


if isnan(maxval) | isnan(minval) | ...
    minval>maxval 
    %errordlg('You must enter a numeric value','Bad Input','modal')
    set(handles.max_edit,'String',strrep(sprintf('%7.3g',handles.oldmax),'',''));
    return
end
%if minval<minimum, minval=minimum; ...
    set(handles.min_edit,'String',strrep(sprintf('%7.3g',minval),'',''));
%end
%if maxval>maximum, maxval=maximum; ... 
    set(handles.max_edit,'String',strrep(sprintf('%7.3g',maxval),'',''));
%end

handles.oldmax=maxval;
handles.oldmin=minval;

set(handles.figure1,'CurrentAxes',handles.axes1);
state=get(handles.dreiD_toggle,'value');

caxis([minval maxval]./(max(handles.data(:))-min(handles.data(:))));
%guidata(h,handles)
%scale_image(handles,0)
update_image(handles,state)
handles.axes5=sp_cbar(handles);
guidata(h,handles)





% --- Executes on button press in pushbutton_flipLR.
function pushbutton_flipLR_Callback(h, eventdata, handles, varargin)
% hObject    handle to pushbutton_flipLR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% We need no reset
set(handles.ToggleReset,'Value',0);

handles.flip = mod(handles.flip + [0 1 0],2);
handles.data = flipdim(handles.data,2);

state=get(handles.dreiD_toggle,'Value');
update_image(handles,state)
handles.axes5=sp_cbar(handles);
guidata(h,handles)


% --- Executes on button press in pushbutton_flipUD.
function pushbutton_flipUD_Callback(h, eventdata, handles, varargin)
% hObject    handle to pushbutton_flipUD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% We need no reset
set(handles.ToggleReset,'Value',0);

handles.flip = mod(handles.flip + [1 0 0],2);
handles.data = flipdim(handles.data,1);

state=get(handles.dreiD_toggle,'Value');
update_image(handles,state)
handles.axes5=sp_cbar(handles);
guidata(h,handles)



% --- Executes on button press in pushbutton_flipBF.
function pushbutton_flipBF_Callback(h, eventdata, handles, varargin)
% hObject    handle to pushbutton_flipBF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.flip = mod(handles.flip + [0 0 1],2);
handles.data = flipdim(handles.data,3);

state=get(handles.dreiD_toggle,'Value');
update_image(handles,state)
handles.axes5=sp_cbar(handles);
guidata(h,handles)

 


function scale_image(handles,state)

maxvalstr=get(handles.max_edit,'String');maxval=str2double(maxvalstr);
minvalstr=get(handles.min_edit,'String');minval=str2double(minvalstr);
map=get(handles.figure1,'colormap');
if state==0, 
ima=handles.data(:,:,1,handles.z);
%updateimage(handles,state)
imghandle=tv(ima,3,[minval maxval],1);
else
ima=squeeze(handles.data(:,:,1,:)); 
imghandle=tv(ima,3,[minval maxval],ceil(sqrt(handles.datasize(3))));
end
set(handles.axes1,'dataaspectratio',handles.aspect)
[img,flag] = getimage(imghandle);
%axis on
%set(gca,'box','on','xtick',[],'ytick',[],'xcolor','b','ycolor','b');
%sp_contrastimage(img,[minval maxval])
sp_contrastimage(handles)

% mgram: Ensure phase view uses a cyclic colormap and wraps endpoints
if get(handles.phase_radio, 'Value') == 1
    colormap(handles.axes1, get_cmp('vikO', 256));
    caxis(handles.axes1, [-pi pi]);   % lock to full phase range
else
    % restore default for non-phase modes
    colormap(handles.axes1, get(0, 'DefaultFigureColormap'));
end
% end mgram

%contrastimage(gcf)
%axis i
%set(gca,'box','on','xtick',[],'ytick',[],'xcolor','b','ycolor','b');
%colormap(map)
%if ishandle(get(handles.axes5,'Children')), colorbar(handles.axes5,'peer',handles.axes1),end

% --------------------------------------------------------------------



% --------------------------------------------------------------------


% --------------------------------------------------------------------
function view_menu_Callback(h, eventdata, handles)
% hObject    handle to view_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function square_submenu_Callback(h, eventdata, handles)
% hObject    handle to square_submenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
state=get(handles.dreiD_toggle,'value');
a=handles.datasize(1)./handles.datasize(2);
handles.aspect=[1 a 1]
guidata(h,handles)
update_image(handles,state);


% --------------------------------------------------------------------
function fourtoone_submenu_Callback(h, eventdata, handles)
% hObject    handle to fourtoone_submenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
state=get(handles.dreiD_toggle,'value');
handles.aspect=[4 1 1];
guidata(h,handles)
update_image(handles,state);

% --------------------------------------------------------------------
function threetoone_submenu_Callback(h, eventdata, handles)
% hObject    handle to threetoone_submenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
state=get(handles.dreiD_toggle,'value');
handles.aspect=[3 1 1];
guidata(h,handles)
update_image(handles,state);

% --------------------------------------------------------------------
function twotoone_submenu_Callback(h, eventdata, handles)
% hObject    handle to twotoone_submenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
state=get(handles.dreiD_toggle,'value');
handles.aspect=[2 1 1];
guidata(h,handles);
update_image(handles,state);


% --------------------------------------------------------------------
function onetoone_submenu_Callback(h, eventdata, handles)
% hObject    handle to onetoone_submenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
state=get(handles.dreiD_toggle,'value');
handles.aspect=[1 1 1];
guidata(h,handles)
update_image(handles,state);

% --------------------------------------------------------------------
function onetotwo_submenu_Callback(h, eventdata, handles)
% hObject    handle to onetotwo_submenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
state=get(handles.dreiD_toggle,'value');
handles.aspect=[1 2 1];
guidata(h,handles)
update_image(handles,state);

% --------------------------------------------------------------------
function onetothree_submenu_Callback(h, eventdata, handles)
% hObject    handle to onetothree_submenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
state=get(handles.dreiD_toggle,'value');
handles.aspect=[1 3 1];
guidata(h,handles)
update_image(handles,state);

% --------------------------------------------------------------------
function onetofour_submenu_Callback(h, eventdata, handles)
% hObject    handle to onetofour_submenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
state=get(handles.dreiD_toggle,'value');
handles.aspect=[1 4 1];
guidata(h,handles)
update_image(handles,state);


% 
% function sp_contrastimage(handles)
% 
% %
% %   Function to provide dynamic adjustment of image contrast
% %
% %   Once the image is displayed, you can hold down a mouse button
% %   and move the mouse to change the image contrast.
% %
% %   Needs helper function updatecontrast.m
% %
% %   Mark Griswold 
% 
% if nargin==1, figHandle=figHandle;else figHandle=gcf;end
% 
% 
% images = findobj(figHandle, 'type', 'image');
% if isempty(images) ,return,end
% im=getimage(images(1));
% if isa(im,'uint8'),return,end
% 
% [ny,nx]=size(im);
% contstruct.clim=[min(im(:)) max(im(:))];
% 
% 
% if isempty(get(findobj(gca,'tag','contdata')))
% 	rmax0 = max(max(im));
% 	rmin0 = min(min(im));
%     %rmin0=0;
% 	contstruct.windowval=0.5;%*(rmax0-rmin0);
% 	contstruct.centerval=contstruct.windowval/2;
%     contstruct.cx=0.5;
%     contstruct.wx=0.5;
% 	contstruct.startpoint=1;
%     contstruct.lastpoint=[ny/2 nx/2 1; ny/2 nx/2 0];
% 	contstruct.scaley=1./(5*ny);
%     contstruct.scalex=1./(5*nx);%(rmax0-rmin0)./100;
% 	axtag=get(gca,'tag');
%     contsruct.ax=gca;
% 	%cla
% 	%imagesc(im);
% 	%axis('square');
%     set(gca,'tag',axtag);
% 	contstruct.axis=get(gca,'children');
%     %contstruct.clim=get(gca,'clim');
%  
% 	t1=text(10,10,' ','tag','contdata','userdata',contstruct,'color','w','fontsize',8);
%     %bdf=get(image,'buttondownfcn')
% 	set(contstruct.axis,'buttondownfcn',['sp_updatecontrast(' '''start''' ')']);
% 	%set(gca,'xgrid','off','xtick',[],'ytick',[]);
%     set(gcf,'doublebuffer','on');
% else
% 	axtag=get(gca,'tag');
% 	contstruct=get(findobj(gca,'tag','contdata'),'userdata');
% 	[n,m]=size(contstruct.axis);
% 	[yy,xx]=size(im);
%     %set(contstruct.axis,'buttondownfcn',['updateimage3d(' '''start''' ')']);
% 	set(contstruct.axis(n),'cdata',im);
% 	kids=get(gca,'children');
% 	set(kids,'visible','off');
% 	set(contstruct.axis(n),'visible','on');
% 	set(gca,'ylim',[0.5 yy+0.5],'xlim',[0.5 xx+0.5]);
%         set(gca,'tag',axtag);
% %	set(findobj(gca,'tag','contdata'),'userdata',contstruct);
% 	set(gca,'xtick',[],'ytick',[]);
% end





% --- Executes on button press in pushbutton_Rotate.
function pushbutton_Rotate_Callback(h, eventdata, handles, varargin)
% hObject    handle to pushbutton_Rotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% We need no reset
set(handles.ToggleReset,'Value',0);

perm = [2 1 3];
handles.data = permute(handles.data,perm);
handles.flip = handles.flip(perm);
handles.perm = handles.perm(perm);
handles.datasize = handles.datasize(perm);
handles.FOV = handles.FOV(perm);


handles.yLim = get(handles.axes1,'Xlim');
handles.xLim = get(handles.axes1,'Ylim');
set(handles.axes1,'Xlim',handles.xLim)
set(handles.axes1,'Ylim',handles.yLim)


state=get(handles.dreiD_toggle,'Value');
update_image(handles,state)
handles.axes5=sp_cbar(handles);
guidata(h,handles)


% --- Executes on button press in pushbutton_Permute.
function pushbutton_Permute_Callback(h, eventdata, handles, varargin)
% hObject    handle to pushbutton_Permute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%perm = handles.perm


% we need no reset
set(handles.ToggleReset,'Value',0);



perm = [3 1 2];
datasize_old = handles.datasize;
z_old       = handles.z;
handles.data = permute(handles.data,perm);
handles.FOV = handles.FOV(perm);
handles.datasize = handles.datasize(perm);
handles.perm = handles.perm(perm);
handles.flip = handles.flip(perm);


handles.z = min(max(ceil(((z_old-1)/(datasize_old(3)-1))*handles.datasize(3)),1),handles.datasize(3));

val=(handles.z-1)/(handles.datasize(3)-1);

set(handles.imgnr_slider,'SliderStep', ...
    [1./(handles.datasize(3)-1) 1./(handles.datasize(3)-1)])

set(handles.imgnr_slider,'Value',val)
imgnrstr=strrep(sprintf('%d',handles.z),' ','');
set(handles.imgnr_edit,'String',imgnrstr)
set(handles.ZoomToggle,'value',0)

state=get(handles.dreiD_toggle,'Value');
update_image(handles,state)
handles.axes5=sp_cbar(handles);
guidata(h,handles)


% --------------------------------------------------------------------
function FOV_submenu_Callback(h, eventdata, handles, varargin)
% hObject    handle to FOV_submenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt={['Enter FOV UD: matrixsize = ' num2str(handles.datasize(1))],...
        ['Enter FOV LR: matrixsize = ' num2str(handles.datasize(2))],...
        ['Enter FOV BF: matrixsize = ' num2str(handles.datasize(3))]};
name='FOV Options';
numlines=1;

defaultanswer={num2str(handles.FOV(1)),num2str(handles.FOV(2)),num2str(handles.FOV(3))};
answer = inputdlg(prompt,name,numlines,defaultanswer);

handles.FOV = [str2double(answer{1}) str2double(answer{2}) str2double(answer{3})];

% no reset
set(handles.ToggleReset,'Value',0);

state=get(handles.dreiD_toggle,'Value');
update_image(handles,state)
guidata(h,handles)
 




% --------------------------------------------------------------------

function update_image(handles,state)

% if get(handles.abs_radio,'value') == 1
%     data = abs(handles.data);
%     1
% end
% 
% if get(handles.real_radio,'value') == 1
%     data = real(handles.data);
%     2
% end
% 
% if get(handles.imag_radio,'value') == 1
%     data = imag(handles.data);
%     3
% end
% 
% if get(handles.phase_radio,'value') == 1
%     data = angle(handles.data);
%     4
% end
% 

axes(handles.axes1)

clim = caxis;
    
reset = get(handles.ToggleReset,'Value');
zoom  = get(handles.ZoomToggle,'Value');

if zoom == 1
    handles.xLim = get(handles.axes1,'xlim');
    handles.yLim = get(handles.axes1,'ylim');
end


if state==0,
    ima=handles.data(:,:,handles.z);
    imghandle=tv(ima,3,[],1);
else
    ima=handles.data;
    imghandle=tv(ima,3,[],ceil(sqrt(handles.datasize(3))));
end

if reset == 1
    clim = [handles.minval handles.maxval];
end

caxis(clim);
set(handles.axes1,'dataaspectratio',handles.FOV./handles.datasize)
%if reset == 0

if zoom == 1
 set(handles.axes1,'Xlim',handles.xLim);
 set(handles.axes1,'ylim',handles.yLim);

end

   %caxis(clim);
   %minval=min(abs(handles.data(:))); minvalstr=strrep(sprintf('%7.3g',minval),'','');
   %maxval=max(abs(handles.data(:))); maxvalstr=strrep(sprintf('%7.3g',maxval),'','');
   %handles.minval=minval;
   %handles.maxval=maxval;
   %set(handles.min_edit,'String',minvalstr)
   %set(handles.max_edit,'String',maxvalstr)
 %  set(handles.ToggleReset,'Value',0);
%end


sp_contrastimage(handles)

% mgram: Ensure phase view uses a cyclic colormap and wraps endpoints
if get(handles.phase_radio, 'Value') == 1
    colormap(handles.axes1, get_cmp('vikO', 256));
    caxis(handles.axes1, [-pi pi]);   % lock to full phase range
else
    % restore default for non-phase modes
    colormap(handles.axes1, get(0, 'DefaultFigureColormap'));
end
% end mgram

%sp_contrastimage(handles)

%set(handles.axes1,'dataaspectratio',handles.aspect)


% --------------------------------------------------------------------


% --- Executes on button press in ToggleReset.
function ToggleReset_Callback(h, eventdata, handles, varargin)
% hObject    handle to ToggleReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ToggleReset
guidata(h, handles);


% --- Executes on button press in ZoomToggle.
function ZoomToggle_Callback(h, eventdata, handles, varargin)
% hObject    handle to ZoomToggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of ZoomToggle
state = get(h,'Value');
zoomout(handles,state)
handles.xLim = get(handles.axes1,'Xlim');
handles.yLim = get(handles.axes1,'ylim');    
%handles.zoom = 1;
guidata(h,handles)
