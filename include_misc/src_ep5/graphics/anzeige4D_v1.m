function CSItoolfig=anzeige4D(daten,limits)
	
% Tool zur anzeige von CSI daten
% und allgemeinen Daten bis maximal 4 Dimensionen
% Akzeptiert als daten Matrix (spekt,Bild1,Bild2,slice)
% 	1.Dim spekt: Spektralerichtung - Richtung die per Mauscursor geplottet wird
% 	2.Dim Bild1 Uebersichtsbild
% 	3.Dim Bild2 Uebersichtsbild
% 	4.Dim slice durchschaltbar ueber Slider
%
% Akzeptiert zudem Struktur:
% 	Feld .xwerte: Vektor(1,n) mit Werten ueber die die Werte aus 1.Dim dargestellt werden sollen
% 	Feld .data_sorted:
% 		zulaessige Datentypen darin:
% 			-Matrix: Sortierung identisch zu vorher (spekt,Bild1,Bild2,slice)
% 			-oder Cell{1,slice}(spekt,Bild1,Bild2)
%			
%
% Written by TBL modfied by VS, modified by TBL, Last Modification February 26, 2009

	% warning off all
	%## Testen in welcher Form die anzuzeigenden Daten uebergeben werden, und uebertragen derselben
	%## in die fuer das Tool gewaehlte Form (Funktion celldat2matdat befindet sich am Ende dieser Datei)

%## Kontrolle ob Daten zur Anzeige uebergeben wurden
	if nargin==0
		%[dateiname,pfad]=uigetfile('*.mat')
		%if dateiname~=0 && strcmpi(dateiname((end-3):end),'.mat')
		%daten=uiimport;
		%if isstruct(daten)
		%	daten=struct2cell(daten);
		%	daten=daten{1};
		%else
		%	return
		%end
		error('Keine anzuzeigenden Daten �bergeben')
	end

%## 'Umformatierung' des Datensatzes auf die vom Tool verwendete Struktur, noetig da verschiedenes als input akzeptiert werden soll
		if ( isstruct(daten)~=1 )
			if ( iscell(daten) )
				datw=celldat2matdat(daten);
			else
				datw=daten;
			end
			clear daten;
			daten.data_sorted=datw;
			clear datw;
			
		elseif ( isfield(daten,'data_sorted') )
			if iscell(daten.data_sorted)
				daten.data_sorted=celldat2matdat(daten.data_sorted);
			end
		
		elseif ( isfield(daten,'data') )
			if ( iscell(daten.data) )
				daten.data_sorted=celldat2matdat(daten.data);
			else
				daten.data_sorted=daten.data;
			end
			daten=rmfield(daten,'data');
		end
		
%## Schreiben einiger zusaetzlich benoetigten Basisgroessen im Datensatz
		[dimdat(1),dimdat(2),dimdat(3),slicenumber]=size(daten.data_sorted);
		spekk0 = abs(permute(sum(sum(daten.data_sorted,3),2),[1 4 2 3]));
		if isreal(daten.data_sorted)
			realwertigeDaten=1;
			abswerte = daten.data_sorted;
		else
			realwertigeDaten=0;
			abswerte = abs(daten.data_sorted);
		end
		
		if isfield(daten,'xwerte')
			xwerte = daten.xwerte;
		else
			xwerte = 1:dimdat(1);
		end
		
%#### erzeugen des Datenfensters mit den Feldern
		
		xclock=clock;
		clockstringfit = [num2str(xclock(4)),'h ',num2str(xclock(5)),'min ',num2str(round(xclock(6))),'s.'];
		
	% Werte initialieren, so dass effektiv ein eintritt in myCallback unterbunden wird ehe die Initialisierung des Skripts beendet ist.
		x11 = 0;	x12 = 0;
		y11 = 0;	y12 = 0;
		x21 = 0;	x22 = 0;
		y21 = 0;	y22 = 0;
		
	fenstersizeskal=[960 690 960 690];		
		
	%## Code um abzufangen das Fenster groesser als Bildschirm ist
	monsize=get(0,'ScreenSize');
	if monsize(1,3)<=(fenstersizeskal(1,3)+20)||monsize(1,4)<=(fenstersizeskal(1,4)+10+30)
		fenstersizeskal_corr=fenstersizeskal.*min([(monsize(1,3)-20)/(fenstersizeskal(1,3)),(monsize(1,4)-10-30)/(fenstersizeskal(1,4))]);
	else
		fenstersizeskal_corr=fenstersizeskal;
	end
		
		CSItoolfig= figure('Name',['CSItooltool opened at: ', clockstringfit],'Position',[10 10 fenstersizeskal_corr(1:2)],'BackingStore','on','WindowButtonDownFcn',@gedrueckt,'Resize','on','BusyAction','cancel','Interruptible','on');
		colormapmenu=uimenu(CSItoolfig,'Label','Colormap');
			uimenu(colormapmenu,'Label','jet','Callback','colormap jet');
			uimenu(colormapmenu,'Label','hot','Callback','colormap hot');
			uimenu(colormapmenu,'Label','hulk','Callback','colormap hulk');
			uimenu(colormapmenu,'Label','ice','Callback','colormap ice');
			uimenu(colormapmenu,'Label','gray','Callback','colormap gray');
			uimenu(colormapmenu,'Label','bone','Callback','colormap bone');
			uimenu(colormapmenu,'Label','copper','Callback','colormap copper');
			uimenu(colormapmenu,'Label','pink','Callback','colormap pink');
			uimenu(colormapmenu,'Label','summer','Callback','colormap summer');
			uimenu(colormapmenu,'Label','hsv','Callback','colormap hsv');
			uimenu(colormapmenu,'Label','cool','Callback','colormap cool');
			uimenu(colormapmenu,'Label','spring','Callback','colormap spring');
			uimenu(colormapmenu,'Label','autumn','Callback','colormap autumn');
			uimenu(colormapmenu,'Label','winter','Callback','colormap winter');
		
		set(gcf,'Color',[0.3,0.7,0.9]);
		%## RecursionLimit begrenzt zahl der nested Funktion aufrufe; Standartwert laut hilfe: large value
		%## teilweise aber seltsame Abbrueche => Test durch setzen auf eine 'sehr grossen Wert'
		set(0,'RecursionLimit',99999999999999999999999999999999999999999999999999999999999);
		
	%## Nur zur optischen Gruppierung der Felder fuer Summenbild
		sumpanel=uipanel('Parent',CSItoolfig,'Title',{'Summenbild'},'Units','normalized','Position',[490 463 90 217]./fenstersizeskal);
		
	%## Initalisieren der Achsen - Positionen fuer die darzustellenden Daten
		axesueb =	axes('Units','normalized','Position',[40 20 400 400]./fenstersizeskal); %Uebersichtsbild
		axesk0 = axes('Units','normalized','Position',[40 540 400 133]./fenstersizeskal);	%Spektrum im K-Raum-Zentrum
		axessum =	axes('Units','normalized','Position',[580 353 320 320]./fenstersizeskal);	%'Integriertes' Uebersichtsbild
		axesspek = axes('Units','normalized','Position',[540 70 400 266]./fenstersizeskal);	%Spektrum im aktuellen Punkt
		
	%## Initialisieren der benoetigten Slider
		if slicenumber>1
			stepslice = 1 / (slicenumber-1);
			sliyafit = uicontrol(CSItoolfig,'Style', 'slider','Units','normalized','Position',[20 20 20 400]./fenstersizeskal,'SliderStep',[stepslice	stepslice],'Callback',@sliderafit);
		else
			stepslice = 1;
			sliyafit = uicontrol(CSItoolfig,'Style', 'slider','Units','normalized','Position',[20 20 20 400]./fenstersizeskal,'SliderStep',[stepslice	stepslice],'Visible','off','Callback',@sliderafit);
		end
		
		if dimdat(1)>1
			spektstep = 1 / (dimdat(1)-1);
			sliybfit = uicontrol(CSItoolfig,'Style', 'slider','Units','normalized','Position',[20 520 440 20]./fenstersizeskal,'SliderStep',[spektstep	spektstep],'Callback',@sliderbfit);
		else
			spektstep = 1;
			sliybfit = uicontrol(CSItoolfig,'Style', 'slider','Units','normalized','Position',[20 520 440 20]./fenstersizeskal,'SliderStep',[spektstep	spektstep],'Visible','off','Callback',@sliderbfit);
		end
		
	%## Felder fuer Koordinaten des Cursors im Uebersichtsbild sowie fuer Slice und Spektrum
		ediuebk1 = uicontrol(CSItoolfig,'Style', 'edit','Units','normalized','Position',[45 480 90 20]./fenstersizeskal,'Callback',@myediuebk12); %slice wert
			uicontrol(CSItoolfig,'Style', 'text','Units','normalized','Position',[45 500 90 10]./fenstersizeskal,'String','Schicht');
		ediuebk2 = uicontrol(CSItoolfig,'Style', 'edit','Units','normalized','Position',[145 480 90 20]./fenstersizeskal,'Callback',@myediuebk12); %spekt wert
			uicontrol(CSItoolfig,'Style', 'text','Units','normalized','Position',[145 500 90 10]./fenstersizeskal,'String','Spekt- Nr');
		ediuebk3 = uicontrol(CSItoolfig,'Style', 'edit','Units','normalized','Position',[245 480 90 20]./fenstersizeskal); %x	wert
			uicontrol(CSItoolfig,'Style', 'text','Units','normalized','Position',[245 500 90 10]./fenstersizeskal,'String','X-Wert');
		ediuebk4 = uicontrol(CSItoolfig,'Style', 'edit','Units','normalized','Position',[345 480 90 20]./fenstersizeskal); %y wert
			uicontrol(CSItoolfig,'Style', 'text','Units','normalized','Position',[345 500 90 10]./fenstersizeskal,'String','Y-Wert');
		
	%## Felder Skalierung des Uebersichtsbildes sowie der Intensitaet im aktuellen Punkt
		ediuebglob = uicontrol(CSItoolfig,'Style','checkbox','Units','normalized','Position',[45 430 90 15]./fenstersizeskal,'String','Skal global','Callback',@myediuebskal);
		ediuebresk = uicontrol(CSItoolfig,'Style','checkbox','Units','normalized','Position',[45 445 90 15]./fenstersizeskal,'String','Reset Skal','Callback',@myediuebskal);
		ediuebmin = uicontrol(CSItoolfig,'Style', 'edit','Units','normalized','Position',[145 430 90 20]./fenstersizeskal,'Callback',@myediuebskal);
			uicontrol(CSItoolfig,'Style','text','Units','normalized','Position',[145 450 90 10]./fenstersizeskal,'String','Skal - min');
		ediuebmax = uicontrol(CSItoolfig,'Style', 'edit','Units','normalized','Position',[245 430 90 20]./fenstersizeskal,'Callback',@myediuebskal);
			uicontrol(CSItoolfig,'Style','text','Units','normalized','Position',[245 450 90 10]./fenstersizeskal,'String','Skal - max');
		ediuebval = uicontrol(CSItoolfig,'Style', 'edit','Units','normalized','Position',[345 430 90 20]./fenstersizeskal);
			uicontrol(CSItoolfig,'Style','text','Units','normalized','Position',[345 450 90 10]./fenstersizeskal,'String','aktueller Wert');
		
	%## Felder zur Skalierung der lokalen Spektrendarstellung
		edispeklok = uicontrol(CSItoolfig,'Style','checkbox','Units','normalized','Position',[540 20 50 15]./fenstersizeskal,'String','lokal','Callback',@myedispekskal);
		edispekmax = uicontrol(CSItoolfig,'Style','checkbox','Units','normalized','Position',[540 35 50 15]./fenstersizeskal,'String','max','Callback',@myedispekskal);
		edispekymin = uicontrol(CSItoolfig,'Style', 'edit','Units','normalized','Position',[590 20 85 20]./fenstersizeskal,'Callback',@myedispekskal); %y wert
			uicontrol(CSItoolfig,'Style', 'text','Units','normalized','Position',[590 40 85 10]./fenstersizeskal,'String','Y-min');
		edispekymax = uicontrol(CSItoolfig,'Style', 'edit','Units','normalized','Position',[675 20 85 20]./fenstersizeskal,'Callback',@myedispekskal); %spekt wert
			uicontrol(CSItoolfig,'Style', 'text','Units','normalized','Position',[675 40 85 10]./fenstersizeskal,'String','Y-max');

		edispekxmin = uicontrol(CSItoolfig,'Style', 'edit','Units','normalized','Position',[770 20 85 20]./fenstersizeskal); %x	wert
			uicontrol(CSItoolfig,'Style', 'text','Units','normalized','Position',[770 40 85 10]./fenstersizeskal,'String','X-min');
		edispekxmax = uicontrol(CSItoolfig,'Style', 'edit','Units','normalized','Position',[855 20 85 20]./fenstersizeskal); %slice wert
			uicontrol(CSItoolfig,'Style', 'text','Units','normalized','Position',[855 40 85 10]./fenstersizeskal,'String','X-max');
		
	%## Felder zur Angabe ueber welche Spektralenpunkte summiert das Integrationsbild erstellt werden soll
		edisumskalmax = uicontrol(CSItoolfig,'Style', 'edit','Units','normalized','Position',[500 473 70 20]./fenstersizeskal,'Callback',@myintsum); %slice wert
			uicontrol(CSItoolfig,'Style', 'text','Units','normalized','Position',[500 493 70 10]./fenstersizeskal,'String','Skal-max');
		edisumskalmin = uicontrol(CSItoolfig,'Style', 'edit','Units','normalized','Position',[500 513 70 20]./fenstersizeskal,'Callback',@myintsum); %slice wert
			uicontrol(CSItoolfig,'Style', 'text','Units','normalized','Position',[500 533 70 10]./fenstersizeskal,'String','Skal-min');
		edisavebut = uicontrol(CSItoolfig,'Style', 'pushbutton','Units','normalized','Position',[500 553 70 30]./fenstersizeskal,'String','SAVE-intsum','Callback',@mysaveintsum);
		edisum1 = uicontrol(CSItoolfig,'Style', 'edit','Units','normalized','Position',[500 593 70 20]./fenstersizeskal,'Callback',@myintsum); %slice wert
			uicontrol(CSItoolfig,'Style', 'text','Units','normalized','Position',[500 613 70 10]./fenstersizeskal,'String','Spekt-Pkts');
		edisumglob = uicontrol(CSItoolfig,'Style','checkbox','Units','normalized','Position',[500 633 70 15]./fenstersizeskal,'String','for all','Callback',@myintsum);
		edisumabs = uicontrol(CSItoolfig,'Style','checkbox','Units','normalized','Position',[500 648 70 15]./fenstersizeskal,'String','abs','Callback',@myintsum);
		
		
	%## Feld zur Aenderung der Dimensionsreihenfolge
		edidimorder = uicontrol(CSItoolfig,'Style', 'edit','Units','normalized','Position',[495 420 70 20]./fenstersizeskal,'Callback',@mypermute);
			uicontrol(CSItoolfig,'Style', 'text','Units','normalized','Position',[495 440 70 10]./fenstersizeskal,'String','Order');
			
	%## Feld zur Bestimmung bei komplexen daten was dargestellt werden soll
		if 	realwertigeDaten==1;
			editypwahl = uicontrol(CSItoolfig,'Style', 'popupmenu','Units','normalized','Position',[495 380 70 20]./fenstersizeskal,'String',{'Real'});
				uicontrol(CSItoolfig,'Style', 'text','Units','normalized','Position',[495 400 70 10]./fenstersizeskal,'String','Data-type');
		else
			editypwahl = uicontrol(CSItoolfig,'Style', 'popupmenu','Units','normalized','Position',[495 380 70 20]./fenstersizeskal,'String',{'Abs','Real','Imag','Phase'},'Callback',@mydatentyptausch);
				uicontrol(CSItoolfig,'Style', 'text','Units','normalized','Position',[495 400 70 10]./fenstersizeskal,'String','Data-type');
		end

		
	%### Initialisieren der Skalierungslisten 
		center = ceil( (dimdat(1)+1)/2 );
		
		skalueb = zeros(2,dimdat(1),slicenumber); % Werte zur Skalierung des Uebersichtsbildes
		intsumskal = zeros(1,2); % Werte zur Skalierung des Selektiven Summen-Uebersichtbildes
		
		maxspekty = repmat([min(abswerte(:)),max(abswerte(:))]',[1,dimdat(2) dimdat(3) slicenumber]);
		maxspekty(2,:,:,:) = maxspekty(2,:,:,:)+((maxspekty(2,:,:,:)-maxspekty(1,:,:,:))==0).*0.00000000000000000001;
		
		intsumskal(1) = min(min(min(abswerte(center,:,:,:),[],4),[],3),[],2);
		intsumskal(2) = max(max(max(abswerte(center,:,:,:),[],4),[],3),[],2);
		
	% Falls in einem Bereich alles gleich, sicherstellen, dass max > min fuer Skalierung
		intsumskal(2) = intsumskal(2)+(intsumskal(1)==intsumskal(2));
		
		skalueb(1,:,:) = squeeze(min(min(abswerte(:,:,:,:),[],3),[],2));
		skalueb(2,:,:) = squeeze(max(max(abswerte(:,:,:,:),[],3),[],2));
		infindex=find(isfinite(skalueb(:))==0);
		skalueb(infindex)=1;
		
	% Falls in einem Bereich alles gleich, sicherstellen, dass max > min fuer Skalierung
		skalueb(2,:,:) = skalueb(2,:,:)+(skalueb(1,:,:)==skalueb(2,:,:));
		
	%### Schreiben der Werte fuer die erste Schicht und Spektrum in die Felder
	% Werte setzen in Hauptfunktion => 'global' innerhalb aller nested Funktionen im Funktionskoerper 
		
		slice=1;
		spekt=1;
		xdata=1;
		ydata=1;
		
	%Erzeugen des Uebersichtsbildes und setzen zugehoeriger Felder
		imghandle=image(permute(abswerte(spekt,:,:,slice),[2 3 4 1]),'CDataMapping','scaled','Parent',axesueb);
		set(axesueb,'CLim',skalueb(:,spekt,slice)')
		axis (axesueb,'off')
		
		set(ediuebk1,'String',num2str(slice))
		set(ediuebk2,'String',num2str(spekt))
		set(ediuebk3,'String','') % Feld fuer xdata
		set(ediuebk4,'String','')	% Feld fuer ydata
		set(ediuebmin,'String',num2str(skalueb(1,spekt,slice)));
		set(ediuebmax,'String',num2str(skalueb(2,spekt,slice)));
        set(get(axesueb,'Parent'),'CurrentAxes',axesueb)
		colorbar('Position',[ 440  20 10 400]./fenstersizeskal);
		
	% Erzeugen des Uebersichtsspektrums (aus k-Raumzentrum) und setzen zugehoeriger Felder
		plot(axesk0,xwerte,spekk0 (:,slice))
		line([xwerte(spekt),xwerte(spekt)],[0,max(spekk0 (:,slice))],'Color','red','Parent',axesk0)
		axis (axesk0,'tight')
		
		% set(edispeklok,'Value',1)
		
		if dimdat(1)==1
			set(edispekxmax,'String',num2str(xwerte(1)+0.1))
		else
			set(edispekxmax,'String',num2str(xwerte(dimdat(1))))
		end
		
		set(edispekxmin,'String',num2str(xwerte(1)))
		set(edispekymin,'String',num2str(maxspekty(1,ydata,xdata,slice)))
		set(edispekymax,'String',num2str(maxspekty(2,ydata,xdata,slice)))
		
	% Erzeugen des Summen Uebersichtsbildes fuer 'SpektrumMittelpunkt' und setzen zugehoeriger Felder
		intmatrix=[]; %Initialisieren der groesse welche die aktuellen summenwerte beinhaltet
		set(edisumglob,'Value',1)
		spektsumlist=repmat({center},[slicenumber,1]);
		spektsumstring=repmat({num2str(center)},[slicenumber,1]);
		
		set(edisumskalmin,'String',intsumskal(1))
		set(edisumskalmax,'String',intsumskal(2))	
		set(edisum1,'String',num2str(spektsumstring{slice}))
		
		myintsumplot
		
	% Schreiben des dimorderfeldes
		if (slicenumber~=1)
			order2=[1 2 3 4];
			set(edidimorder,'String','1 2 3 4')%,'BackgroundColor',[0,0,1],'ForegroundColor',[1,1,1])
		else
			order2=[1 2 3];
			set(edidimorder,'String','1 2 3')%,'BackgroundColor',[0,0,1],'ForegroundColor',[1,1,1])
		end
		
		handles.pointer = 'arrow';
		
	if nargin==2
	set(ediuebmin,'String',num2str(limits(1)))
	set(ediuebmax,'String',num2str(limits(2)))
	set(ediuebglob,'Value',1);
	myediuebskal
	end
	
	
	%## Setzen einiger Werte, die in myCallback verwendet werden
		
	%Skalierungsfaktor fuer Uebersichtsbild um aus der Position den darunterliegenden Pixel zu ermitteln
		xdatafactor1 = (dimdat(3)/(400/fenstersizeskal(1)));
		ydatafactor1 = (dimdat(2)/(400/fenstersizeskal(2)));
		
	%Skalierungsfaktor fuer summiertes Uebersichtsbild um aus der Position den darunterliegenden Pixel zu ermitteln
		xdatafactor2 = (dimdat(3)/(320/fenstersizeskal(1)));
		ydatafactor2 = (dimdat(2)/(320/fenstersizeskal(2)));
	
	%Bereich indem das Uebersichtsbild dargestellt wird - [40 20 400 400]./fenstersizeskal
		x11 = 40/fenstersizeskal(1);	x12 = 440/fenstersizeskal(1);
		y11 = 20/fenstersizeskal(2);	y12 = 420/fenstersizeskal(2); 
		
	%Bereich indem das summierte Uebersichtsbild dargestellt wird - [580 353 320 320]./fenstersizeskal
		x21 = 580/fenstersizeskal(1);	x22 = 900/fenstersizeskal(1);
		y21 = 353/fenstersizeskal(2);	y22 = 673/fenstersizeskal(2); 
	
	
%################## eingebettete Funktion myCallback - wird bei gedrueckter Maustaste bei jeder Mausbewegung aufgerufen
	
function gedrueckt(src,eventdata)
		set(src,'WindowButtonMotionFcn',@myCallback)
		set(src,'WindowButtonUpFcn',@geloest)
		myCallback
	end
	
function geloest(src,eventdata)
		set(src,'Pointer','arrow')
		set(src,'WindowButtonMotionFcn','')
		set(src,'WindowButtonUpFcn','')
	end
	
function myCallback(src,eventdata)
		
		% Abfragen der Position von der Maus und dem CSItoolfenster und umrechnen in Relative Einheiten
		xy=get(0,'PointerLocation');
		Position=get(CSItoolfig,'Position');
		xy1=(xy-Position(1:2))./Position(3:4);
		
	% Abfrage ob Maus im Bereich des Uebersichtsbildes
		if xy1(1) > x11 && xy1(1) < x12 && xy1(2) > y11 && xy1(2) < y12
			
			set(CSItoolfig,'Pointer','crosshair')
			
	% uebersetzten der xy Werte in Pixel Variablen
	% pruefen ob Maus auf neuen Pixel gerutscht ist - nur wenn ja aktualisieren
			
			xdatn=ceil(xdatafactor1 * (xy1(1)- x11));
			ydatn=dimdat(2)-ceil(ydatafactor1 * (xy1(2)-y11))+1;
			if ((xdata ~= xdatn)||(ydata ~= ydatn))
				
	% Fragt den gewuenschten Auschschnitsbereich ab - test ob Werte Zahlen sind
				xmin1 = str2double(get(edispekxmin,'String'));
				xmax1 = str2double(get(edispekxmax,'String'));
				if isnan([xmin1 xmax1])==0
					xmin=xmin1;
					xmax=xmax1;
				end
				
				xdata = xdatn;
				ydata = ydatn;
				
				set(ediuebk3,'String',num2str(xdata)) % x-wert
				set(ediuebk4,'String',num2str(ydata)) % y-wert
				set(ediuebval,'String',num2str(abswerte(spekt,ydata,xdata,slice))) % Wert in Voxel
				set(edispekymax,'String',num2str(maxspekty(2,ydata,xdata,slice)))
				set(edispekymin,'String',num2str(maxspekty(1,ydata,xdata,slice)))
				set(edispekxmax,'String',xmax)
				set(edispekxmin,'String',xmin)
				
				plot(axesspek,xwerte,abswerte(:,ydata,xdata,slice),'blue');
				line([xwerte(spekt),xwerte(spekt)],[maxspekty(:,ydata,xdata,slice)'],'Color','red','Parent',axesspek);
				axis(axesspek,[xmin xmax (maxspekty(:,ydata,xdata,slice)')])
				set(axesspek,'XGrid','on','YGrid','on')
				
				%drawnow
				
			end
			
	% Abfrage ob Maus im Bereich des summierten Uebersichtsbildes
		elseif xy1(1) > x21 && xy1(1) < x22 && xy1(2) > y21 && xy1(2) < y22
			
			set(CSItoolfig,'Pointer','crosshair')
			
			% uebersetzten der xy Werte in Pixel Variablen
			% pruefen ob Maus auf neuen Pixel gerutscht ist - nur wenn ja aktualisieren
			
			xdatn=ceil(xdatafactor2 * (xy1(1)- x21));
			ydatn=dimdat(2)-ceil(ydatafactor2 * (xy1(2)-y21))+1;
			if ((xdata ~= xdatn)||(ydata ~= ydatn))
				
				% Fragt den gewuenschten Auschschnitsbereich ab - test ob Werte Zahlen sind
				xmin1 = str2double(get(edispekxmin,'String'));
				xmax1 = str2double(get(edispekxmax,'String'));
				if isnan([xmin1 xmax1])==0
					xmin=xmin1;
					xmax=xmax1;
				end
			
				xdata = xdatn;
				ydata = ydatn;
				
				set(ediuebk3,'String',num2str(xdata)) % x-wert
				set(ediuebk4,'String',num2str(ydata)) % y-wert
				set(ediuebval,'String',num2str(intmatrix(ydata,xdata))) % Wert in Voxel
				set(edispekymax,'String',num2str(maxspekty(2,ydata,xdata,slice)))
				set(edispekymin,'String',num2str(maxspekty(1,ydata,xdata,slice)))
				set(edispekxmax,'String',xmax)
				set(edispekxmin,'String',xmin)
												
				plot(axesspek,xwerte,abswerte(:,ydata,xdata,slice),'blue');
				line([xwerte(spekt),xwerte(spekt)],[maxspekty(:,ydata,xdata,slice)'],'Color','red','Parent',axesspek);
				axis(axesspek,[xmin xmax (maxspekty(:,ydata,xdata,slice)')])
				set(axesspek,'XGrid','on','YGrid','on')
				
				%drawnow
				
			end
		else
			set(CSItoolfig,'Pointer','arrow')	
		end
		
	end
	
	
%############## Eingebettete Funktion sliderafit - aktualisiert nach umschalten der Slice das Uebersichtsbild und passt Wert slice an.
	
function	sliderafit(src,eventdata)
		
		slice=round((get(sliyafit,'Value'))/stepslice+1);
		
		set(ediuebmin,'String',num2str(skalueb(1,spekt,slice)));
		set(ediuebmax,'String',num2str(skalueb(2,spekt,slice)));
		
		imghandle=image(permute(abswerte(spekt,:,:,slice),[2 3 4 1]),'CDataMapping','scaled','Parent',axesueb);
		set(axesueb,'CLim',skalueb(:,spekt,slice)');
		axis (axesueb,'off')
		set(get(axesueb,'Parent'),'CurrentAxes',axesueb)
        colorbar('Position',[ 440  20 10 400]./fenstersizeskal);
		
		plot(axesk0,xwerte,spekk0 (:,slice))
		line([xwerte(spekt),xwerte(spekt)],[0,max(spekk0(:,slice))],'Color','red','Parent',axesk0)
		axis (axesk0,'tight')
		
		set(ediuebk1,'String',num2str(slice));
		
	% aktualisieren des Summiertenuebersichtsbildes - aufruf der zugehoerigen Funktion
		set(edisum1,'String',spektsumstring{slice})
		myintsumplot
		
	end
	
	
%############## Eingebettete Funktion sliderbfit - Aendert Position im Spektrum und passt Uebersichtsbild an.
	
function	sliderbfit(src,eventdata)
		
		spekt=round((get(sliybfit,'Value'))/spektstep+1);
		
		set(ediuebmin,'String',num2str(skalueb(1,spekt,slice)));
		set(ediuebmax,'String',num2str(skalueb(2,spekt,slice)));
		
		imghandle=image(permute(abswerte(spekt,:,:,slice),[2 3 4 1]),'CDataMapping','scaled','Parent',axesueb);
		set(axesueb,'CLim',skalueb(:,spekt,slice)');
		axis (axesueb,'off')
		set(get(axesueb,'Parent'),'CurrentAxes',axesueb)
        colorbar('Position',[ 440  20 10 400]./fenstersizeskal);
		
		plot(axesk0,xwerte,spekk0(:,slice))
		line([xwerte(spekt),xwerte(spekt)],[0,max(spekk0(:,slice))],'Color','red','Parent',axesk0)
		axis (axesk0,'tight')
		
		set(ediuebk2,'String',num2str(spekt))
			
	end
	
	
%############## Eingebettete Funktion myediuebk12 - ermoeglicht aenderung der Spektposition und Sliceposition durch Zahleingabe
	
function	myediuebk12(src,eventdata)
		
		spekt1=str2double(get(ediuebk2,'String'));
		slice1=str2double(get(ediuebk1,'String'));
		
		if (ismember(spekt1,1:dimdat(1)) && ismember(slice1,1:slicenumber))
			
			spekt=spekt1;
			newslice=slice-slice1; % Dient dazu um feststellen zu koennen ob myintsum aufgerufen werden muss
			slice=slice1;
			set(ediuebmin,'String',num2str(skalueb(1,spekt,slice)));
			set(ediuebmax,'String',num2str(skalueb(2,spekt,slice)));
			
			imghandle=image(permute(abswerte(spekt,:,:,slice),[2 3 4 1]),'CDataMapping','scaled','Parent',axesueb);
			set(axesueb,'CLim',skalueb(:,spekt,slice)');
			axis (axesueb,'off')
			set(get(axesueb,'Parent'),'CurrentAxes',axesueb)
            colorbar('Position',[ 440  20 10 400]./fenstersizeskal);
			
			plot(axesk0,xwerte,spekk0(:,slice))
			line([xwerte(spekt),xwerte(spekt)],[0,max(spekk0(:,slice))],'Color','red','Parent',axesk0)
			axis (axesk0,'tight')
			
%			set(ediuebk2,'String',[num2str(spekt)])
			set(sliybfit,'Value',(spekt-1)*spektstep);
			set(sliyafit,'Value',(slice-1)*stepslice);
			
			% aktualisieren des Summiertenuebersichtsbildes - aufruf der zugehoerigen Funktion falls andere Slice
			if (newslice~=0)
				myintsumplot
			end
			
		else
			set(ediuebk2,'String',num2str(spekt))
			set(ediuebk1,'String',num2str(slice))
		end
		
	end
	
	
%############## Eingebettete Funktion myskal - aendern der Skalierungswerte das Uebersichtsbild und reset auf Startwerte
	
function	myediuebskal(src,eventdata)
		
		uebglob=get(ediuebglob,'Value');
		uebresk=get(ediuebresk,'Value');
		
		if uebresk==1
			skalueb(1,:,:)=squeeze(min(min(abswerte(:,:,:,:),[],3),[],2));
			skalueb(2,:,:)=squeeze(max(max(abswerte(:,:,:,:),[],3),[],2));
			skalueb(2,:,:)=skalueb(2,:,:)+(skalueb(1,:,:)==skalueb(2,:,:));
			set(ediuebresk,'Value',0)
			set(ediuebglob,'Value',0)
		else
			skmin=str2double(get(ediuebmin,'String'));
			skmax=str2double(get(ediuebmax,'String'));
						
			if (isnan([skmin,skmax])==0)
				if skmin < skmax
					if uebglob==1;
						for c1=1:slicenumber
							for c2=1:dimdat(1)
								skalueb(1:2,c2,c1)=[skmin,skmax];
							end
						end
					else
						skalueb(:,spekt,slice)=[skmin,skmax];
					end
				end
			end
		end
		
		imghandle=image(permute(abswerte(spekt,:,:,slice),[2 3 4 1]),'CDataMapping','scaled','Parent',axesueb);
		set(axesueb,'CLim',skalueb(:,spekt,slice)');
		axis (axesueb,'off')
		set(get(axesueb,'Parent'),'CurrentAxes',axesueb)
        colorbar('Position',[ 440  20 10 400]./fenstersizeskal);
		
		set(ediuebmin,'String',num2str(skalueb(1,spekt,slice)))
		set(ediuebmax,'String',num2str(skalueb(2,spekt,slice)))
		
	end
	
	
%############# Eingebettete Funktion myedispekskal - setzen der y-Skalierung im Spektrum	
	%Definieren dieser beiden werte um ueberpruefen zu koennen welcher button zuletzt gedrueckt wurde
	speklok=0;
	spekmax=0;
	
function myedispekskal(src,eventdata)
		
	%Sicherstellen, das immer nur ein Button aktiv ist, da sie sich gegenseitig ausschliesen
		speklok1 = get(edispeklok,'Value');
		spekmax1 = get(edispekmax,'Value');
		maxs1=str2double(get(edispekymax,'String'));
		mins1=str2double(get(edispekymin,'String'));
		
		% Pruefen ob ymin und ymax Zahlenwerte enthalten, wenn nein rueckschreiben der urspruenglichen Werte
		% au�erdem ist im letzteren Fall klar, dass deswegen diese unterfunktion ausgefuehrt wird und somit ein neu schreiben der Listen unnoetig ist
		if ( (isnan(mins1)==0) && (isnan(maxs1)==0) && (maxs1 > mins1) )
			maxs=maxs1;
			mins=mins1;
			nurskal=0;
		else
			maxs=maxspekty(2,ydata,xdata,slice);
			mins=maxspekty(1,ydata,xdata,slice);
			nurskal=1;
		end
		
	% Falls beide Felder markiert sind pruefen welches zuletzt aktiviert wurde und das andere zurueck auf 0 setzen
		if (speklok1==1 && spekmax1==1)
			if speklok==0
				speklok1 = 1;
				spekmax1 = 0;
				nurskal=0;
				
			else
				speklok1 = 0;
				spekmax1 = 1;
				nurskal=0;
			end
	% und ausserdem testen ob ein checkbox wert sich veraendert hat, wenn nein dann Listen potentiell nur lokal veraendern
		elseif ( (speklok==speklok1) && (spekmax == spekmax1))
			nurskal=1;
		end
		
		speklok=speklok1;
		spekmax=spekmax1;
		
	% Schreiben der skalierungswerte
		if speklok==1
			if nurskal==0
				maxspekty = [(min(abswerte(:,:,:,:),[],1));(max(abswerte(:,:,:,:),[],1))];
				maxspekty(2,:,:,:) = maxspekty(2,:,:,:)+((maxspekty(2,:,:,:)-maxspekty(1,:,:,:))==0).*0.00000000000000000001;
			else
				maxspekty(:,ydata,xdata,slice)=[mins;maxs];
			end
			
		elseif ( spekmax==1 && nurskal==0)
			maxspekty = repmat([min(abswerte(:)),max(abswerte(:))]',[1,dimdat(2) dimdat(3) slicenumber]);
			maxspekty(2,:,:,:) = maxspekty(2,:,:,:)+((maxspekty(2,:,:,:)-maxspekty(1,:,:,:))==0).*0.00000000000000000001;
			
		else
			maxspekty = repmat([mins;maxs],[1,dimdat(2) dimdat(3) slicenumber]);
		end
		
		set(edispeklok,'Value',speklok)
		set(edispekmax,'Value',spekmax)
		set(edispekymax,'String',num2str(maxspekty(2,ydata,xdata,slice)))
		set(edispekymin,'String',num2str(maxspekty(1,ydata,xdata,slice)))
		
	end
	
	
%############## Eingebettete Funktion myintsum - daten fuer Integrationsummenbild aufbereiten und bilderstellung aufrufen

function	myintsum(src,eventdata)
		
		pktallslice = get(edisumglob,'Value');
		sumskal=zeros(1,2);
		intsum1=get(edisum1,'String');
		intmask=zeros(dimdat(1),1);
		
	% Zerteilen der Feldwerte an Kommapositionen und schreiben der Werte in eine Liste
		intsum=textscan(intsum1,'%s','Delimiter',',');
		intsum=intsum{1};
		anzahl=length(intsum);
		
	% Ersetzen vorhandener bis zeichen ('-') durch ':'
		for c1=1:anzahl
			intsum{c1}=strrep(intsum{c1},'-',':');
		end
		
	% Uebertragen der Werte in eine Auswahlliste
		intwerte=[];
		for c1=1:anzahl
			if sum(eval(['ismember(',intsum{c1},',intwerte)']))~=0
				errordlg('einzelne Spektrale Punkte mehrfach - Jeder wird nur einmal genommen','Fehler','modal');
			end
			eval(['intwerte=union(',intsum{c1},',intwerte);']);
		end
		
	% Schreiben der Werte fuer diese oder alle Schichten
		if (pktallslice==1)
			spektsumlist=repmat({intwerte},[slicenumber,1]);
			spektsumstring=repmat({intsum1},[slicenumber,1]);
		else
			spektsumlist(slice)={intwerte};
			spektsumstring(slice)={intsum1};
		end
		
	% Abfragen der Skalierungswerte und pruefen ob sie Zahlen enthalten
		sumskal(1)=str2double(get(edisumskalmin,'String'));
		sumskal(2)=str2double(get(edisumskalmax,'String'));
		
		if (isnan(sumskal)==0)
			intsumskal(1:2)=sumskal(:);
		end
		myintsumplot
	end

	
%############## Eingebettete Funktion myintsumplot - Integrationsummenbild erstellen

function	myintsumplot(src,eventdata)
		
		sumabsdat = get(edisumabs,'Value');
		
	% Addition der abs oder rohdaten und Normierung durch Anzahl der Punkte
		if sumabsdat==1;
			% Versuch unn�tige Berechnungen zu vermeiden: Falls komplex und datentyp Abs dann summiere die abswerte, sonst arbeite auf orginal daten
			if ( realwertigeDaten==0 && get(editypwahl,'Value')==1)
				intmatrix=permute(sum(abswerte(spektsumlist{slice},:,:,slice),1)./length(spektsumlist{slice}),[2 3 4 1]);
			else
				intmatrix=permute(sum(abs(daten.data_sorted(spektsumlist{slice},:,:,slice)),1)./length(spektsumlist{slice}),[2 3 4 1]);
			end
		else
			intmatrix=permute(sum(daten.data_sorted(spektsumlist{slice},:,:,slice),1)./length(spektsumlist{slice}),[2 3 4 1]);
		end
		
		set(edisumskalmin,'String',intsumskal(1))
		set(edisumskalmax,'String',intsumskal(2))
		
		if isreal(intmatrix)
			image(intmatrix,'CDataMapping','scaled','Parent',axessum);
			set(axessum,'CLim',intsumskal(1,:))
		else
			image(abs(intmatrix),'CDataMapping','scaled','Parent',axessum);
			set(axessum,'CLim',intsumskal(1,:))
		end
		axis (axessum,'off')
	end
	
%############## Eingebettete Funktion mysaveintsum - Integrationsummenbild erstellen und skalieren
	
function	mysaveintsum(src,eventdata)
		sumabsdat = get(edisumabs,'Value');
		
		answer = mysaveinterface;
		
		if ( isempty(answer)~=1 )
			if sumabsdat==1;
				for c1=1:slicenumber
					intmatrixs(c1,:,:)=sum(abs(daten.data_sorted(spektsumlist{c1},:,:,c1)),1)./length(spektsumlist{c1});
				end
			else
				for c1=1:slicenumber
					intmatrixs(c1,:,:)=sum(daten.data_sorted(spektsumlist{c1},:,:,c1),1)./length(spektsumlist{c1});
				end
			end
			
			savefct(intmatrixs,spektsumstring,answer{1},answer{2},answer{3});
		end
		
	end

	
%############## Eingebettete Funktion mydatentyptausch - Einstellen des zu zeigenden Datentypes bei komplexen Daten: ABS, REAL, ...
% wird nur aufgerufen falls uebergebene Daten komplexwertig

function	mydatentyptausch(src,eventdata)
		Feldnummer = get(editypwahl,'Value');
		Auswahlakt = get(editypwahl,'String');
		
	% Vermeiden, dass waehrend die Daten ueberarbeitet werden myCallback ausgefuehrt wird
		x11 = 0;	x12 = 0;
		y11 = 0;	y12 = 0;
		x21 = 0;	x22 = 0;
		y21 = 0;	y22 = 0;
		
		switch Auswahlakt{Feldnummer}
			case 'Abs'
				abswerte = abs(daten.data_sorted);
			case 'Real'
				abswerte = real(daten.data_sorted);
			case 'Imag'
				abswerte = imag(daten.data_sorted);
			case 'Phase'
				abswerte = angle(daten.data_sorted);
		end
			
	%### Skalierungslisten setzen - abhaengig von den gedrueckten Skalierungsbuttons
		uebglob1 = get(ediuebglob,'Value');
		speklok1 = get(edispeklok,'Value');
		
		if uebglob1==0
			skalueb=zeros(2,dimdat(1),slicenumber); % Werte zur Skalierung des Uebersichtsbildes
			skalueb(1,:,:)=squeeze(min(min(abswerte(:,:,:,:),[],3),[],2));
			skalueb(2,:,:)=squeeze(max(max(abswerte(:,:,:,:),[],3),[],2));
			
			% Falls in einem Bereich alles gleich, sicherstellen, dass max > min fuer Skalierung
			skalueb(2,:,:)=skalueb(2,:,:)+(skalueb(1,:,:)==skalueb(2,:,:));
		else
			skalueb=repmat(skalueb(:,1,1),[1 dimdat(1) slicenumber]);
		end
		
		maxspekty = repmat([min(abswerte(:)),max(abswerte(:))]',[1,dimdat(2) dimdat(3) slicenumber]);
		maxspekty(2,:,:,:) = maxspekty(2,:,:,:)+((maxspekty(2,:,:,:)-maxspekty(1,:,:,:))==0).*0.00000000000000000001;
			
		
	%### Schreiben der Werte fuer die akt. Schicht und Spektrum in die Felder
	% Erzeugen des Uebersichtsbildes und setzen zugehoeriger Felder
		imghandle=image(permute(abswerte(spekt,:,:,slice),[2 3 4 1]),'CDataMapping','scaled','Parent',axesueb);
		set(axesueb,'CLim',skalueb(:,spekt,slice)');
		axis (axesueb,'off')
		set(get(axesueb,'Parent'),'CurrentAxes',axesueb)
        colorbar('Position',[ 440  20 10 400]./fenstersizeskal);
		
		set(ediuebmin,'String',num2str(skalueb(1,spekt,slice)));
		set(ediuebmax,'String',num2str(skalueb(2,spekt,slice)));
		
	% Erzeugen des Uebersichtsspektrums (aus k-Raumzentrum) und setzen zugehoeriger Felder
		plot(axesk0,xwerte,spekk0 (:,slice))
		line([xwerte(spekt),xwerte(spekt)],[0,max(spekk0 (:,slice))],'Color','red','Parent',axesk0)
		axis (axesk0,'tight')
		
		set(edispekymin,'String',num2str(maxspekty(1,ydata,xdata,slice)))
		set(edispekymax,'String',num2str(maxspekty(2,ydata,xdata,slice)))
		
	% Erzeugen des Summen Uebersichtsbildes fuer 'SpektrumMittelpunkt' und setzen zugehoeriger Felder
		myintsum
		
	% Ruecksetzen der aktiven Fensterbereiche auf 'Betriebswerte'
		x11 =  40/fenstersizeskal(1);	x12 = 440/fenstersizeskal(1);
		y11 =  20/fenstersizeskal(2);	y12 = 420/fenstersizeskal(2); 
		x21 = 580/fenstersizeskal(1);	x22 = 900/fenstersizeskal(1);
		y21 = 353/fenstersizeskal(2);	y22 = 673/fenstersizeskal(2); 
			
	end
	
%############## Eingebettete Funktion mypermute - Permutieren der Dimensionen des Datensatzes
	
function	mypermute(src,eventdata)
	
	% Abfragen noetiger Werte
		order1=get(edidimorder,'String');
		uebglob1 = get(ediuebglob,'Value');
		speklok1 = get(edispeklok,'Value');
		
		order=eval(['[',order1,']']);
		
	% Erzeugen eines effektiven Umsortierungsvektors, da daten nicht mehr gemaess urspruenglicher
		% Sortierung vorliegen muessen
		for c1=1:length(order)
			[k,order3(c1)]=max(order(c1)==order2);
		end
		
	% Vermeiden, dass waehrend die Daten umsortiert werden myCallback ausgefuehrt wird
		x11 = 0;	x12 = 0;
		y11 = 0;	y12 = 0;
		x21 = 0;	x22 = 0;
		y21 = 0;	y22 = 0;
		
	% Abfrage ob der fall gegeben ist, bei dem spekk0 nicht neu berechnet werden muss
		daten.data_sorted=permute(daten.data_sorted,order3);
		if ( order(1)==order2(1) && length(order)==4 && order(4)==order2(4))
		else
			spekk0 = abs(permute(sum(sum(daten.data_sorted,3),2),[1 4 2 3]));
		end
		
		order2=order;
		
	% Umsortieren des Datensatzes
		abswerte = permute(abswerte,order3);
		[dimdat(1),dimdat(2),dimdat(3),slicenumber]=size(abswerte);
		
		if ( isfield(daten,'xwerte') && order(1)==1 )
				xwerte = daten.xwerte;
		else
				xwerte = 1:dimdat(1);
		end
		
		slice=1;
		spekt=1;
		xdata=1;
		ydata=1;
		
	% Skalierungsfaktoren um die darunterliegenden Pixel zu ermitteln
		xdatafactor1 = (dimdat(3)/(400/940));
		ydatafactor1 = (dimdat(2)/(400/690));
		xdatafactor2 = (dimdat(3)/(320/940));
		ydatafactor2 = (dimdat(2)/(320/690));
		
	% Slider an/aus
		if slicenumber>1
			stepslice = 1 / (slicenumber-1);
			set(sliyafit,'Visible','on','Value',(slice-1)*stepslice,'SliderStep',[stepslice	stepslice])
		else
			stepslice = 1;
			set(sliyafit,'Visible','off')
		end
		
		if dimdat(1)>1
			spektstep = 1 / (dimdat(1)-1);
			set(sliybfit,'Visible','on','Value',(spekt-1)*spektstep,'SliderStep',[spektstep	spektstep])
		else
			spektstep = 1;
			set(sliybfit,'Visible','off')
		end	
		
	%### Skalierungslisten setzen - abhaengig von den gedrueckten Skalierungsbuttons
		center=ceil((dimdat(1)+1)/2);
		
		if uebglob1==0
			skalueb=zeros(2,dimdat(1),slicenumber); % Werte zur Skalierung des Uebersichtsbildes
			skalueb(1,:,:)=squeeze(min(min(abswerte(:,:,:,:),[],3),[],2));
			skalueb(2,:,:)=squeeze(max(max(abswerte(:,:,:,:),[],3),[],2));
			
			% Falls in einem Bereich alles gleich, sicherstellen, dass max > min fuer Skalierung
			skalueb(2,:,:)=skalueb(2,:,:)+(skalueb(1,:,:)==skalueb(2,:,:));
		else
			skalueb=repmat(skalueb(:,1,1),[1 dimdat(1) slicenumber]);
		end
		
		if speklok1==1
			maxspekty = [(min(abswerte(:,:,:,:),[],1));(max(abswerte(:,:,:,:),[],1))];
			maxspekty(2,:,:,:) = maxspekty(2,:,:,:)+((maxspekty(2,:,:,:)-maxspekty(1,:,:,:))==0).*0.00000000000000000001;
		else
			maxspekty=repmat(maxspekty(:,1,1,1),[1 dimdat(2) dimdat(3) slicenumber]);
		end
			
		
	%### Schreiben der Werte fuer die erste Schicht und Spektrum in die Felder
	% Erzeugen des Uebersichtsbildes und setzen zugehoeriger Felder
		imghandle=image(permute(abswerte(spekt,:,:,slice),[2 3 4 1]),'CDataMapping','scaled','Parent',axesueb);
		set(axesueb,'CLim',skalueb(:,spekt,slice)');
		axis (axesueb,'off')
		set(get(axesueb,'Parent'),'CurrentAxes',axesueb)
        colorbar('Position',[ 440  20 10 400]./fenstersizeskal);
		
		set(ediuebk1,'String',num2str(slice))
		set(ediuebk2,'String',num2str(spekt))
		set(ediuebk3,'String',num2str(xdata))
		set(ediuebk4,'String',num2str(ydata))
		set(ediuebmin,'String',num2str(skalueb(1,spekt,slice)));
		set(ediuebmax,'String',num2str(skalueb(2,spekt,slice)));
		
	% Erzeugen des Uebersichtsspektrums (aus k-Raumzentrum) und setzen zugehoeriger Felder
		plot(axesk0,xwerte,spekk0 (:,slice))
		line([xwerte(spekt),xwerte(spekt)],[0,max(spekk0 (:,slice))],'Color','red','Parent',axesk0)
		axis (axesk0,'tight')
		
		if dimdat(1)==1
			set(edispekxmax,'String',num2str(xwerte(1)+0.1))
		else
			set(edispekxmax,'String',num2str(xwerte(dimdat(1))))
		end
		set(edispekxmin,'String',num2str(xwerte(1)))
		set(edispekymin,'String',num2str(maxspekty(1,ydata,xdata,slice)))
		set(edispekymax,'String',num2str(maxspekty(2,ydata,xdata,slice)))
		
	% Erzeugen des Summen Uebersichtsbildes fuer 'SpektrumMittelpunkt' und setzen zugehoeriger Felder
		spektsumlist=repmat({center},[slicenumber,1]);
		spektsumstring=repmat({num2str(center)},[slicenumber,1]);
		set(edisum1,'String',num2str(spektsumstring{1}))
		myintsum
		
	% Ruecksetzen der aktiven Fensterbereiche auf 'Betriebswerte'
		x11 =  40/fenstersizeskal(1);	x12 = 440/fenstersizeskal(1);
		y11 =  20/fenstersizeskal(2);	y12 = 420/fenstersizeskal(2); 
		x21 = 580/fenstersizeskal(1);	x22 = 900/fenstersizeskal(1);
		y21 = 353/fenstersizeskal(2);	y22 = 673/fenstersizeskal(2); 
	
	end
	
end % Schliessen der aeusseren funktion

%### Alle folgenden Funktion sind auserhalb des Hauptfunktionskoerpers! - d.h. sie kennen nur die ihnen uebergebenen Variablen
%        d.h. keine Gefahr von Namenskollisionen!!!


% Funktion mysavedinterface erzeugt Fenster zur Datenabfrage beim speichern
function answer=mysaveinterface
	
	b=get(0,'ScreenSize');
	pos=round(b(3:4)/2)-[100 80];
	mydialog=dialog('Name','Speichern','Position',[pos 200 160],'CloseRequestFcn',@mycancel);
	
	fenstersizeskalsave=[200 160 200 160];
	edivarname = uicontrol(mydialog,'Style', 'edit','Units','normalized','Position',[ 15 115 170 20]./fenstersizeskalsave,'String','intsum','BackgroundColor','w','HorizontalAlignment','left');
			uicontrol(mydialog,'Style', 'text','Units','normalized','Position',[ 30 135 90 10]./fenstersizeskalsave,'String','Variablenname','HorizontalAlignment','left');
	edifilename = uicontrol(mydialog,'Style', 'edit','Units','normalized','Position',[ 15 75 170 20]./fenstersizeskalsave,'String','csidata.mat','BackgroundColor','w','HorizontalAlignment','left');
			uicontrol(mydialog,'Style', 'text','Units','normalized','Position',[ 30 95 90 10]./fenstersizeskalsave,'String','Dateiname','HorizontalAlignment','left');
	ediappendsave = uicontrol(mydialog,'Style', 'checkbox','Units','normalized','Position',[ 30 45 70 20]./fenstersizeskalsave,'String','Append?');
	
	edisavebutton = uicontrol(mydialog,'Style', 'pushbutton','Units','normalized','Position',[ 105 15 70 20]./fenstersizeskalsave,'String','Save','Callback',@mysave);
	edicancelbutton = uicontrol(mydialog,'Style', 'pushbutton','Units','normalized','Position',[ 15 15 70 20]./fenstersizeskalsave,'String','Cancel','Callback',@mycancel);
	
	uiwait(mydialog)
	
	function mysave(src,eventdata)
		varname=get(edivarname,'String');
		filename=get(edifilename,'String');
		appendyes=get(ediappendsave,'Value');
		answer={varname,filename,num2str(appendyes)};
		delete (mydialog)
	end

	function mycancel(src,eventdata)
		delete (mydialog)
		answer={};
		
	end

end

% Funktion zum umbenennen und speichern der Bilddaten - ausgelagert zur Vermeidung von Namenskollisionen!!!
function savefct(intmatrix,spektsumstring,variablenname,dateiname,hinzufuegen)
		eval([variablenname,'.data_sorted=intmatrix;'])
		eval([variablenname,'.auswahl=spektsumstring;'])
		hinzufuegen=str2double(hinzufuegen);
		if hinzufuegen==1
			save (dateiname,variablenname,'-append');
		else 
			save (dateiname,variablenname);
		end
end


% Funktion um Daten die Schichtweise in Cellen einsortiert sind in eine 4D-Matrix einzusortieren
function matdaten=celldat2matdat(celldaten)
	slicenumber=length(celldaten);
	dimdat=size(celldaten{1});
	matdaten=permute(reshape([celldaten{:}],dimdat(1),dimdat(2),slicenumber,dimdat(3)),[1,2,4,3]);
	
	% alte version
	%slicenumber=length(celldaten);
	%dimdat=size(celldaten{1});
	%matdaten=zeros(dimdat(1),dimdat(2),dimdat(3),slicenumber);
	%for cx1=1:slicenumber
	%	matdaten(:,:,:,cx1)=celldaten{cx1};
	%end
end