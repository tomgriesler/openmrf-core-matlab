%% DesignPhantom.m
%
% GUI to draw a phantom composed of ellipses, polygons and Bezier curves.
%
% INPUT: * (optional) a background image
%
% OUTPUT: * the phantom structure as used by the other functions in the
%           toolbox.
%
% Notes:
% 1) The origin is at the center of the FOV
% 2) the spatial coordinates (center, width, vertex, controls points)
% are defined in FOV units; i.e., with a FOV of 25 cm the point (1,-1) is
% located at x=25cm and y=-25cm to the origin.
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 10-02-2011 (dd-mm-yyyy)

function phantom = DesignPhantom(background,phantom)
close all;
if  nargin==0
    background = zeros(128);
end
if nargin<2
   phantom.region = [];
end
if max(abs(background(:)))>0
    background = double(background);
    background = abs(background)/max(abs(background(:)));
end
res = size(background);
control = [];
hnodes = [];
hcurve = [];
indRegion = max(0,numel(phantom.region)-1);
str = [];
x = (-0.5:0.1:0.5);
shift = .5*(1+mod(res+1,2));
while numel(str)==0
    tmp = inputdlg('Field of view in your chosen unit (must be coherent with k-space positions):','Please answer',1,{'1'},struct('WindowStyle','normal'));
    if isempty(tmp)
        str = [];
    else
        str = tmp{1};
    end
end
phantom.FOV = str2double(str)*[1,1];
im = RasterizePhantom(phantom,res,1);

message_standard = 'Standard mode: Press ''h'' for help, ''q'': stop';
message_modif = 'Modification mode: Left click to select point to modify, right click to end mode';
scrsz = get(0,'ScreenSize');
fig2 = figure('OuterPosition',[600,scrsz(4),600,600]);axes('DrawMode','fast');imagesc(im);axis image;colormap gray;colorbar;
set(gca,'XTick',(x+0.5)*res(1),'XTickLabel',x*phantom.FOV(1),'YTick',(x+0.5)*res(2),'YTickLabel',x*phantom.FOV(2),'XDir','normal','YDir','reverse');
fig1 = figure('OuterPosition',[0,scrsz(4),600,600]);ah = axes('DrawMode','fast');imagesc(background-im);axis image;colormap jet;colorbar;drawnow;
set(gca,'XTick',(x+0.5)*res(1),'XTickLabel',x*phantom.FOV(1),'YTick',(x+0.5)*res(2),'YTickLabel',x*phantom.FOV(2),'XDir','normal','YDir','reverse');
title(message_standard);
type = {'Ellipse','Polygon','Bezier'};

[Ntype,DrawRegion,weight] =  new_region();

quit = false;
while ~quit % wait for the user to quit GUI
    pause(0.1);
end

if nargout==0    
    ExportPhantom(phantom,'NewPhantom');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions (GUI)

    function [Ntype,DrawRegion,weight] = new_region()
        figure(fig1);
        indRegion = indRegion+1;
        Ntype = menu('Type of region',type);
        weight = UIAskWeight();
        switch Ntype
            case 1 % ELLIPSE
                DrawRegion = @DrawEllipse;
            case 2 % POLYGON
                DrawRegion = @DrawPolygon;
            case 3 % BEZIER
                DrawRegion = @DrawBezier;
        end
        zoom out;
        set(fig1,'WindowButtonDownFcn',@wbdf,'WindowButtonMotionFcn',@wbmf,'WindowKeyPressFcn',@wkpf);
        cp = get(ah,'CurrentPoint');
        pos = [cp(1,1) cp(1,2)]-1;
        [hnodes,hcurve] = DrawRegion(pos);
        
        function val = UIAskWeight()
            str = [];
            while numel(str)==0
                tmp = inputdlg('Weight of next region:','Please answer',1,{'1'},struct('WindowStyle','normal'));
                if isempty(tmp)
                    str = [];
                else
                    str = tmp{1};
                end
            end
            val = str2double(str);
        end
    end

    function region = convert_to_region(control,Ntype,weight)
        region.type = lower(type{Ntype});
        region.weight = weight;
        x = (control(:,2)/res(1)-1/2);% center and scale the axis
        y = (control(:,1)/res(2)-1/2);
%        if ispolycw(x, y)&&(Ntype>1) %ensure that the control points are defined counter-clockwise for polygons and Bezier
%            [x, y] = poly2ccw(x,y);
%        end
        switch Ntype
            case 1 % ELLIPSE
                % 1: the center, 2: a border point, 3: a focal point
                region.center = [x(1),y(1)];
                v = 2*([x(1)-x(3),y(1)-y(3)]);
                region.angle = atan2(v(2),v(1));
                L = sqrt((x(2)-x(3))^2+(y(2)-y(3))^2)+sqrt((x(2)-2*x(1)+x(3))^2+(y(2)-2*y(1)+y(3))^2);
                region.width = 2*[L/2, real(sqrt(L^2-v(:)'*v(:))/2)];
            case 2 % POLYGON
                region.vertex(:,1) = x;
                region.vertex(:,2) = y;
            case 3 % BEZIER
                region.control(:,1) = x;
                region.control(:,2) = y;
                %flag = CheckRegion(region);
        end
    end

    function wbdf(src,~) % When pressing the mouse
        switch get(src,'SelectionType')
            case 'normal' % LEFT CLICK
                cp = get(ah,'CurrentPoint');
                pos = [cp(1,1) cp(1,2)]-1;
                control = [control; pos];
                [hnodes,hcurve] = DrawRegion(control,hnodes,hcurve);title(message_standard);
            case 'alt' % RIGHT CLICK
                if size(control,1)>2
                    phantom.region{indRegion} = convert_to_region(control,Ntype,weight);
                    phtmp.FOV = phantom.FOV;phtmp.region = cell(1,1);phtmp.region{1} = phantom.region{indRegion};
                    im = im + RasterizePhantom(phtmp,res,1);
                    figure(fig2);imagesc(im);axis image;colorbar;title('current phantom');
                    figure(fig1);imagesc(background-im);axis image;colorbar;title('current phantom');DrawRegion(control);
                    set(src,'WindowButtonDownFcn','','WindowButtonMotionFcn','','WindowKeyPressFcn','');
                    control = [];
                    hnodes = [];
                    hcurve = [];
                    [Ntype,DrawRegion,weight] = new_region();
                else
                    disp('You should validate at least 3 control points before finishing with the element.');
                end
                cp = get(ah,'CurrentPoint');
                pos = [cp(1,1) cp(1,2)]-1;
                [hnodes,hcurve] = DrawRegion([control;pos],hnodes,hcurve);title(message_standard);
            case 'extend' % EXTENDED CLICK
                disp('you performed an ''extended'' click');
            otherwise
                disp('action to be precised');
        end
    end

    function wbmf(~,~,ind) % MOVING THE MOUSE
        cp = get(ah,'CurrentPoint');
        pos = [cp(1,1) cp(1,2)]-1;
        if nargin<3
            [hnodes,hcurve] = DrawRegion([control; pos],hnodes,hcurve);title(message_standard);
            %text(pos(1),pos(2),sprintf('(%0.1f,%0.1f)',pos(1),pos(2)),'VerticalAlignment','bottom');drawnow
        else
            cc = control;
            cc(ind,:) = pos;
            [hnodes,hcurve] = DrawRegion(cc,hnodes,hcurve);title(message_modif);
        end
    end

    function wkpf(src,evnt) % PRESSING A KEY
        %disp(['you pressed key: ' evnt.Key]);
        switch evnt.Key
            case 'm' % MODIFY
                disp('enter modification mode');
                if numel(control)>0
                    [hnodes,hcurve] = DrawRegion(control,hnodes,hcurve);
                    title(message_modif);
                    set(src,'WindowButtonDownFcn',@wbdf_mod);
                    set(src,'WindowButtonMotionFcn','');
                end
            case 'r' % RASTERIZE PHANTOM
                disp('Rendering the phantom');
                phtmp = phantom;
                if size(control,1)>2
                    phtmp.region{indRegion} = convert_to_region(control,Ntype,weight);
                end
                figure(fig2);RasterizePhantom(phtmp,res,1);
                figure(fig1);
            case 'q' % QUIT
                disp('about to finish');
                set(src,'WindowButtonDownFcn','');
                set(src,'WindowButtonMotionFcn','');
                if size(control,1)>2
                    phantom.region{indRegion} = convert_to_region(control,Ntype,weight);
                    phtmp.FOV = phantom.FOV;phtmp.region = cell(1,1);phtmp.region{1} = phantom.region{indRegion};
                    im = im + RasterizePhantom(phtmp,res,1);
                    figure(fig2);imagesc(im);axis image;colorbar;
                    figure(fig1);DrawRegion(control,hnodes,hcurve);
                else
                    warning('PhantonGUI:quit','You did not validate enough control points for the last element.');
                    if ~exist('phantom','var')
                        phantom = [];
                    end
                end
                quit = true;
                close(fig1);
            case 'z' % UNDO
                cp = get(ah,'CurrentPoint');
                pos = [cp(1,1) cp(1,2)]-1;
                control = control(1:end-1,:);
                [hnodes,hcurve] = DrawRegion([control; pos],hnodes,hcurve);
            case 'h'
                disp('GUI for drawing Bezier-based phantoms:');
                disp('- Standard mode:');
                disp('  * Left click:  add control point');
                disp('  * Right click: finish current curve and start new one');
                disp('  * ''q'': finish current curve and stop drawing');
                disp('  * ''m'': enter modification mode');
                disp('  * ''z'': undo last control point');
                disp('  * ''r'': render the current phantom');
                disp('  * ''h'': display this help');
                disp('- Modification mode:');
                disp('  * Left click and hold:  select and move a control point');
                disp('  * Right click: leave modification mode');
            otherwise
                evnt.Key
        end
        
        function wbdf_mod(src,~) % When pressing the mouse in modification mode
            switch get(src,'SelectionType')
                case 'normal' % LEFT CLICK
                    cp = get(ah,'CurrentPoint');
                    pos = [cp(1,1) cp(1,2)]-1;
                    distances = sqrt(abs(control(:,1)-pos(1)).^2+abs(control(:,2)-pos(2)).^2);
                    [~, ind_selected] = min(distances);
                    title('Modification mode: move the point and right click');
                    set(src,'WindowButtonMotionFcn',{@wbmf,ind_selected},'WindowButtonUpFcn',{@wbuf,ind_selected});
                case 'alt' % RIGHT CLICK
                    disp('enter standard mode');
                    set(src,'WindowButtonDownFcn',@wbdf,'WindowButtonMotionFcn',@wbmf,'WindowButtonUpFcn','');
            end
        end
        
        function wbuf(src,~,ind) % When releasing the mouse in modification mode
            cp = get(ah,'CurrentPoint');
            pos = [cp(1,1) cp(1,2)]-1;
            control(ind,:) = pos;
            switch get(src,'SelectionType')
                case 'normal' % LEFT CLICK
                    title('Modification mode: move the point and right click');
                    set(src,'WindowButtonDownFcn',@wbdf_mod);
                    set(src,'WindowButtonMotionFcn','');
                case 'alt' % RIGHT CLICK
                    disp('enter standard mode');
                    set(src,'WindowButtonDownFcn',@wbdf);
                    set(src,'WindowButtonMotionFcn',@wbmf);
                    set(src,'WindowButtonDownFcn',@wbdf,'WindowButtonMotionFcn',@wbmf,'WindowButtonUpFcn','');
            end
        end
    end

    function [hnodes,hcurve] = DrawBezier(curve,hnodes,hcurve)
        %% DrawBezier.m
        % Draws a polygon that corresponds to a fine approximation of the
        % specified quadratic Bezier curve.
        % Input:    * N control points in a Nx2 matrix.
        if nargin<3, hcurve = [];end
        if nargin<2, hnodes = [];end
        curve(:,1) = curve(:,1)+shift(1);% the shift ensures consistency with rasterized phantoms (see: RasterizePhantom.m)
        curve(:,2) = curve(:,2)+shift(2);
        if isempty(hnodes)
            hold on;hnodes = plot(curve(:,1),curve(:,2),'y+-','MarkerSize',7);hold off;
        else
            set(hnodes,'XData',curve(:,1),'YData',curve(:,2));
        end
        N = size(curve,1);
        if N>1
            dt = 1/10;
            t = 0:dt:N;
            p = repmat((0:N-1),length(t),1);
            t = repmat(t',1,N);
            M = beta2(t-p);M(:,end-1:end) = M(:,end-1:end)+beta2(t(:,1:2)-repmat([-2,-1],[size(t,1),1]));
            f = M*curve;
        else
            f = curve;
        end
        if isempty(hcurve)
            hold on;hcurve = plot(f(:,1),f(:,2),'m-','Linewidth',1);hold off;
        else
            set(hcurve,'XData',f(:,1),'YData',f(:,2));
        end
        drawnow;
    end
    function [hnodes,hcurve] = DrawPolygon(curve,hnodes,hcurve)
        %% DrawPolygon.m
        % Draws the specified polygon
        % Input:    * N vertices of the polygon in a Nx2 matrix.
        if nargin<3, hcurve = [];end
        if nargin<2, hnodes = [];end
        curve(:,1) = curve(:,1)+shift(1);% the shift ensures consistency with rasterized phantoms (see: RasterizePhantom.m)
        curve(:,2) = curve(:,2)+shift(2);
        if isempty(hnodes)
            hold on;hnodes = plot(curve(:,1),curve(:,2),'y+','MarkerSize',7);hold off;
        else
            set(hnodes,'XData',curve(:,1),'YData',curve(:,2));
        end
        if isempty(hcurve)
            hold on;hcurve = plot([curve(:,1); curve(1,1)],[curve(:,2); curve(1,2)],'m-','Linewidth',1);hold off;
        else
            set(hcurve,'XData',[curve(:,1); curve(1,1)],'YData',[curve(:,2); curve(1,2)]);
        end
        drawnow;
    end
    function [hnodes,hcurve] = DrawEllipse(M,hnodes,hcurve)
        %% DrawEllipse.m
        % Draws the specified ellipse
        % Input:    * (3,2) matrix that defines three points
        %               (the center, a border point, and a focal point)
        if nargin<3, hcurve = [];end
        if nargin<2, hnodes = [];end
        M(:,1) = M(:,1)+.5*(1+mod(res(1)+1,2));% the shift ensures consistency with rasterized phantoms (see: RasterizePhantom.m)
        M(:,2) = M(:,2)+.5*(1+mod(res(2)+1,2));
        if isempty(hnodes)
            hold on;hnodes = plot(M(:,1),M(:,2),'y+','MarkerSize',7);hold off;
        else
            set(hnodes,'XData',M(:,1),'YData',M(:,2));
        end
        n = min(3,size(M,1));
        switch n
            case 0
                error('not enought points in the input');
            case 1
                f = M(1,:).';
            case 2
                center = M(1,:);
                m = M(2,:);
                R = norm(center-m,2);
                ct = cos((0:2*pi/30:2*pi));
                st = sin((0:2*pi/30:2*pi));
                f = [R*ct+center(1);R*st+center(2)];
            case 3
                center = M(1,:);
                m = M(2,:);
                f1 = M(3,:);
                f2 = 2*center-f1;
                v = f1-f2;
                angle = atan2(v(2),v(1));
                L = norm(f1-m,2)+norm(f2-m,2);
                B = real(sqrt(L^2-v(:)'*v(:))/2);
                A = real(L/2);
                ct = cos((0:2*pi/50:2*pi));
                ct0 = cos(angle);
                st = sin((0:2*pi/50:2*pi));
                st0 = sin(angle);
                f = [A*ct*ct0-B*st*st0+center(1);B*st*ct0+A*ct*st0+center(2)];
        end
        if isempty(hcurve)
            hold on;hcurve = plot([f(1,:), f(1,1)],[f(2,:), f(2,1)],'m-','Linewidth',1);hold off;
        else
            set(hcurve,'XData',[f(1,:), f(1,1)],'YData',[f(2,:), f(2,1)]);
        end
        hold off;
        drawnow;
    end
end