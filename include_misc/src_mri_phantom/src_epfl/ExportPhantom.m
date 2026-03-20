
%% ExportPhantom.m
%
% This function export a phantom into files with the following format:
%           * SVG(Z):   visualize it with the Opera web browser
%           * TEX:      using TIKZ, will produce clean PDF using pdflatex
%           * M:        Calling this script in matlab defines the phantom
%
% INPUTS:   * phantom:  structure that defines the phantom
%           * filename: name of the phantom
%           * Decimals: max. number of decimals to be used in the file
%           * PIXEL_WIDTH (fixe this value to 0 for a resizable SVG)
%           * GAMMA (gamma correction exponent for SVG. A value different than 1 could cause improper display correctly on some viewers)
%
% SEE: DEFINESL, DESIGNPHANTOM
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 26-03-2011 (dd-mm-yyyy)

function ExportPhantom(phantom,filename,Decimals,PIXEL_WIDTH,GAMMA)
%% Default arguments
if nargin<5, GAMMA = 1;end
if nargin<4, PIXEL_WIDTH = 0;end
if nargin<3, Decimals = 4;end
if nargin<2, filename = 'phantom';end
if nargin<1
    phantom.FOV = [1,1];
    phantom.region = cell(4,1);
    phantom.region{1} = struct( 'type', 'ellipse', 'weight', 0.2, 'center', [0 0], 'angle', 0.0000, 'width', [0.85 0.85]);
    phantom.region{2} = struct( 'type', 'polygon', 'weight', 0.4, 'vertex', [0.26 -0.15; 0.06 0.4; -0.4 0.08]);
    phantom.region{3} = struct( 'type', 'bezier', 'weight', 0.3, 'control', [-0.06 -0.27; 0.11 0.04; -0.23 0.31; -0.34 -0.27]);
    phantom.region{4} = struct( 'type', 'polygon', 'weight', -0.2, 'vertex', [-0.14 0.08; -0.15 -0.11; -0.02 0.05]);
    phantom.region{5} = struct( 'type', 'polygon', 'weight', 0.1, 'vertex', [-.2 .1; -.15 .1; -.2 .05]);
end
% Use Bruno Luong's insidepoly.m
% (see http://www.mathworks.com/matlabcentral/fileexchange/27840-2d-polygon-interior-detection)
% if possible or fallback to matlab's inpolygon
if exist('insidepoly_dblengine','file')==3, inpoly = @insidepoly;
else inpoly = @inpolygon;end

%% Headers
svgfile = fopen([filename '.svg'],'w+');
fprintf(svgfile,'<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n');
fprintf(svgfile,'<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n');
if ~(PIXEL_WIDTH==0)
    fprintf(svgfile,'<svg xmlns:svg="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns="http://www.w3.org/2000/svg" version="1.1" viewBox="0 0 %dpx %dpx">\n',PIXEL_WIDTH,PIXEL_WIDTH);
else
    fprintf(svgfile,'<svg xmlns:svg="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns="http://www.w3.org/2000/svg" version="1.1" preserveAspectRatio="xMinYMin meet" viewBox="-.5 -.5 1 1">\n');
end
fprintf(svgfile,'<title>%s</title>\n',filename);
fprintf(svgfile,'<desc>Generated with ExportPhantom.m on %s</desc>\n<desc>Recommanded viewers: recent verion of Opera or Firefox</desc>\n',datestr(clock()));
fprintf(svgfile,'<desc>Contact: Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne</desc>\n');

texfile = fopen([filename '.tex'],'w+');
fprintf(texfile,'\\documentclass{article}\n\\usepackage{tikz}\n\\usepackage[active,tightpage]{preview}\n\\PreviewEnvironment{tikzpicture}\n\\setlength{\\PreviewBorder}{0bp}\n');
%fprintf(texfile,'\\usepackage[pdftex]{hyperref}\n\\hypersetup{pdftitle={%s phantom},pdfcreator={Generated with ExportPhantom.m},pdfproducer={Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne},pdfstartview={FitB}}\n',filename);
fprintf(texfile,'\\begin{document}\n\\begin{tikzpicture}\n');
%fprintf(texfile,'\\usetikzlibrary{calc}\n\\newcommand{\\qbz}[3]{(#1)..controls($(#1)!{2/3}!(#2)$)and($(#2)!{1/3}!(#3)$)..}\n');
headerstr = sprintf('%% %s\n%% Generated with ExportPhantom.m on %s\n%% Contact: Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne\n',upper(filename),datestr(clock()));
fprintf(texfile,'%s',headerstr);

mfile = fopen([filename '.m'],'w+');
name  = filename(regexp(lower(filename),'[a-z]'));
fprintf(mfile,'%s',headerstr);
fprintf(mfile,'%s.FOV = [%.5g,%.5g];\n',name,phantom.FOV(1),phantom.FOV(2));

%% Definitions
fprintf(svgfile,'<defs>\n');
if GAMMA~=1
    fprintf(svgfile,'\t<filter id="Gamma"><feComponentTransfer><feFuncR type="gamma" exponent="%s"/><feFuncG type="gamma" exponent="%s"/><feFuncB type="gamma" exponent="%s"/></feComponentTransfer></filter>\n',strNum(GAMMA),strNum(GAMMA),strNum(GAMMA));
end
N = numel(phantom.region);
indexDigits = floor(log10(N)+1);

% declare the regions
fprintf(mfile,'%s.region = cell(%d,1);\n',name,N);
for idRegion = 1:N
    region = phantom.region{idRegion};
    switch region.type
        case 'ellipse'
            if (region.width(1)==region.width(2))
                fprintf(svgfile,'\t<circle  id="%s" cx="%s" cy="%s" r="%s"/>\n',strRegion(idRegion),strNum(region.center(2)),strNum(region.center(1)),strNum(region.width(1)/2));
                fprintf(texfile,'\\def\\%s{(%s,%s)circle(%s)}\n',strRegion(idRegion),strNum(region.center(2)),strNum(-region.center(1)),strNum(region.width(1)/2));
            else
                if region.angle==0
                    fprintf(svgfile,'\t<ellipse id="%s" rx="%s" ry="%s" transform="translate(%s,%s)"/>\n',...
                        strRegion(idRegion),strNum(region.width(2)/2),strNum(region.width(1)/2),strNum(region.center(2)),strNum(region.center(1)));
                    fprintf(texfile,'\\def\\%s{(%s,%s)ellipse(%s and %s)}\n',strRegion(idRegion),strNum(region.center(2)),strNum(-region.center(1)),strNum(region.width(2)/2),strNum(region.width(1)/2));
                else
                    fprintf(svgfile,'\t<ellipse id="%s" rx="%s" ry="%s" transform="translate(%s,%s) rotate(%s,0,0)"/>\n',...
                        strRegion(idRegion),strNum(region.width(2)/2),strNum(region.width(1)/2),strNum(region.center(2)),strNum(region.center(1)),strNum(-region.angle*180/pi));
                    fprintf(texfile,'\\def\\%s{[shift={(%s,%s)},rotate=%s](0,0)ellipse(%s and %s)}\n',strRegion(idRegion),strNum(region.center(2)),strNum(-region.center(1)),strNum(region.angle*180/pi),strNum(region.width(2)/2),strNum(region.width(1)/2));
                end
            end
            fprintf(mfile,'%s.region{%d} = struct(''type'',''ellipse'',''weight'',%s,''center'',%s,''angle'',%s,''width'',%s);\n',name,idRegion,strNum(region.weight),vertex2m(region.center),strNum(region.angle),vertex2m(region.width));
        case 'polygon'
            fprintf(svgfile,'\t<polygon id="%s" points="%s"/>\n',strRegion(idRegion),vertex2svg(region.vertex));
            fprintf(texfile,'\\def\\%s{%scycle}\n',strRegion(idRegion),vertex2tex(region.vertex,'--'));
            fprintf(mfile,'%s.region{%d} = struct(''type'',''polygon'',''weight'',%s,''vertex'',%s);\n',name,idRegion,strNum(region.weight),vertex2m(region.vertex));
        case 'bezier'
            svgpath = control2svg(region.control);
            fprintf(svgfile,'\t<path    id="%s" d="%s"/>\n',strRegion(idRegion),svgpath);
            fprintf(texfile,'\\def\\%s{%s}\n',strRegion(idRegion),control2tex(region.control));
            % with PGF version>=2.1:
            %fprintf(texfile,'\\def\\%s{[scale=25,yscale=-1] svg "%s"}\n',strRegion(idRegion),svgpath);
            fprintf(mfile,'%s.region{%d} = struct(''type'',''bezier'',''weight'',%s,''control'',%s);\n',name,idRegion,strNum(region.weight),vertex2m(region.control));
        otherwise
            error('unknown region type');
    end
end
% detect structure of the phantom
approx = cell(1,N);
for idRegion = 1:N
    region = phantom.region{idRegion};
    switch region.type
        case 'ellipse'
            A = real(region.width(1)/2);B = real(region.width(2)/2);
            dt = 2*pi/60;
            t = 0:dt:2*pi;
            ct = cos(t);ct0 = cos(region.angle);
            st = sin(t);st0 = sin(region.angle);
            tmp = [A*ct*ct0-B*st*st0+region.center(1);B*st*ct0+A*ct*st0+region.center(2)].';
            s = [0;cumsum(sqrt(sum((tmp(1:end-1,:)-tmp(2:end,:)).^2,2)))].';
            ds = s(end)/60;
            t = interp1(s,t,0:ds:s(end));
            ct = cos(t);ct0 = cos(region.angle);
            st = sin(t);st0 = sin(region.angle);
            approx{idRegion} = [A*ct*ct0-B*st*st0+region.center(1);B*st*ct0+A*ct*st0+region.center(2)].';
        case 'polygon'
            approx{idRegion} = region.vertex([1:end,1],:);
        case 'bezier'
            dt = 1/15;
            nodes = control2node(region.control);
            Nnodes = size(nodes,1);
            t = 0:dt:Nnodes;
            p = repmat((0:Nnodes-1),length(t),1);
            t = repmat(t',1,Nnodes);
            M = beta2(t-p);M(:,end-1:end) = M(:,end-1:end)+beta2(t(:,1:2)-repmat([-2,-1],[size(t,1),1]));
            approx{idRegion} = (M*region.control);
%             tmp = (M*region.control);
%             s = [0;cumsum(sqrt(sum((tmp(1:end-1,:)-tmp(2:end,:)).^2,2)))].';
%             ds = s(end)/10/Nnodes;
%             t = interp1(s,0:dt:Nnodes,0:ds:s(end));
%             p = repmat((0:Nnodes-1),length(t),1);
%             t = repmat(t',1,Nnodes);
%             M = beta2(t-p);M(:,end-1:end) = M(:,end-1:end)+beta2(t(:,1:2)-repmat([-2,-1],[size(t,1),1]));
%             approx{idRegion} = (M*region.control).';
            clear M p t;
    end
end
level = 1;
intersection{level}.elements = (1:N)';
intensity = zeros(1,N);
for id=1:N
    intensity(id) = phantom.region{id}.weight;
end
intersection{level}.intensity = intensity;
intersection{level}.priority = zeros(1,N);
Ni = N;
maxintens = 0;
minintens = 0;
maxpriority = 0;
minpriority = 0;
while Ni>0
    intersection{level+1}.elements = [];
    intersection{level+1}.intensity = [];
    intersection{level+1}.priority = [];
    
    testing_regions = unique(sort(intersection{level}.elements(:)))';
    %fprintf('****** LEVEL %d *******\n',level);
    for idtest = testing_regions
        v_test = approx{idtest};
        for idcur = find(all(intersection{level}.elements~=idtest,2))'
            cur = intersection{level}.elements(idcur,:);
            in = ones(1,size(v_test,1));
            new_intens = intersection{1}.intensity(idtest);
            k=1;
            is_cur_in_test = zeros(1,length(cur));
            is_test_in_cur = zeros(1,length(cur));
            for idregion = cur
                v_cur = approx{idregion};
                ind = find(in);
                in = 0*in;
                in(ind(inpoly(v_test(ind,1),v_test(ind,2),v_cur(:,1),v_cur(:,2)))) = 1;
                is_cur_in_test(k) = numel(find(inpoly(v_cur(:,1),v_cur(:,2),v_test(:,1),v_test(:,2))))==size(v_cur,1);
                is_test_in_cur(k) = numel(find(inpoly(v_test(:,1),v_test(:,2),v_cur(:,1),v_cur(:,2))))==size(v_test,1);
                new_intens = new_intens+phantom.region{idregion}.weight;
                k = k+1;
            end
            if all(is_cur_in_test) % the current element is included is the test region
                %fprintf('reg %d encloses element %d (regions %s)\n',idtest,idcur,sprintf('%d ',cur));
                intersection{1}.priority(idtest) = intersection{1}.priority(idtest)+1;
                maxpriority = max(maxpriority,intersection{1}.priority(idtest));
                %fprintf('    -> increasing priority of reg %d (level 1); new value: %d\n',idtest,intersection{1}.priority(idtest));
            end
            if numel(find(in))>0 %
                if numel(find(in))==size(v_test,1)
                    %fprintf('reg %d is included in element %d (regions %s)\n',idtest,idcur,sprintf('%d ',cur));
                    %intersection{level}.priority(idcur) = intersection{level}.priority(idcur)+1;
                    intersection{1}.priority(idtest) = intersection{1}.priority(idtest)-1;
                    minpriority = min(minpriority,intersection{1}.priority(idtest));
                    intersection{level}.intensity(idtest) = new_intens;%sum(intersection{level}.intensity([idtest,idcur]));
                    %fprintf('   -> decreasing priority of reg %d (level 1); new value: %d\n',idtest,intersection{1}.priority(idtest));
                    %fprintf('   -> element %d (regions %s) gets intensity %.2f\n',idtest,sprintf('%d ',intersection{level}.elements(idtest,:)),new_intens);
                else
                    %fprintf('reg %d intersects element %d (regions %s)\n',idtest,idcur,sprintf('%d ',cur));
                    if ~any(is_test_in_cur)
                        new_elements = unique([intersection{level+1}.elements;sort([cur,idtest])],'rows');
                        if  size(new_elements,1)>size(intersection{level+1}.elements,1)
                            intersection{level+1}.elements = new_elements;
                            %fprintf('   -> creating a new element (regions %s) for level %d\n',sprintf('%d ',sort([cur,idtest])),level+1);
                        end
                    end
                end
            end
        end
    end
    Ni = size(intersection{level+1}.elements,1);
    intersection{level+1}.priority = zeros(1,Ni);
    intersection{level+1}.intensity = zeros(1,Ni);
    for idcur = 1:Ni
        cur = intersection{level+1}.elements(idcur,:);
        new_intens = intersection{1}.intensity(cur(end));
        for idcur2 = cur(1:end-1)
            new_intens = new_intens+phantom.region{idcur2}.weight;
        end
        intersection{level+1}.intensity(idcur) = new_intens;
        %fprintf('   -> new element (regions %s) gets intensity %.2f\n',sprintf('%d ',cur),new_intens);
        intersection{level+1}.priority(idcur) = min(intersection{1}.priority(cur));
        %fprintf('   -> new element (regions %s) gets priority %d\n',sprintf('%d ',cur),intersection{level+1}.priority(idcur));
    end
    maxintens = max([intersection{level}.intensity,intersection{level+1}.intensity,maxintens]);
    minintens = min([intersection{level}.intensity,intersection{level+1}.intensity,minintens]);
    level = level+1;
end
maxintens
for level=2:numel(intersection)-1
    elements=unique(intersection{level}.elements(:,1:end-1));
    for idel=1:size(elements,1)
        el = elements(idel,:);
        switch numel(el)
            case 1
                fprintf(svgfile,'\t<clipPath id="X%s"><use xlink:href="#%s"/></clipPath>\n',strRegion(el(1)),strRegion(el(1)));
            otherwise
                id = [];
                for reg=el(1:end-1)
                    id = sprintf('%sX%s',id,strRegion(reg));
                end
                fprintf(svgfile,'\t<clipPath id="%s" clip-path="url(#%s)"><use xlink:href="#%s"/></clipPath>\n',sprintf('%sX%s',id,strRegion(el(end))),id,strRegion(el(end)));
        end
    end
end
fprintf(svgfile,'</defs>\n');

%% Display elements
gstr = [];
if (PIXEL_WIDTH~=0)
    gstr = sprintf('%s transform="scale(%d,%d) translate(0.5,0.5)"',gstr,PIXEL_WIDTH,PIXEL_WIDTH);
end
if GAMMA~=1
    gstr = sprintf('%s filter="url(#Gamma)"',gstr);
end
fprintf(svgfile,'<g%s>\n\t<rect width="1" height="1" x="-.5" y="-.5"/>\n',gstr);
fprintf(texfile,'\\begin{scope}[scale=128,shift={(.5,.5)}]\n\\color[gray]{0}\\fill(-.5,-.5)rectangle(.5,.5);\n');

for priority = maxpriority:-1:minpriority
    %fprintf('***** PRIORITY %d ******\n',priority);
    for level=1:numel(intersection)-1
        elements=intersection{level}.elements;
        for idel=1:size(elements,1)
            if intersection{level}.priority(idel)==priority
                el = elements(idel,:);
                if numel(el)==1
                    intensity = max(0,100*intersection{level}.intensity(idel)/maxintens);
                    %fprintf('region: %s, with intensity: %.1f\n',sprintf('%d ',el),intensity/100);
                    fprintf(svgfile,'\t<use xlink:href="#%s" fill="rgb(%s%%,%s%%,%s%%)"/>\n',strRegion(el),strNum(intensity),strNum(intensity),strNum(intensity));
                    %fprintf(svgfile,'\t<use xlink:href="#%s" fill="hsl(0,0%%,%s%%)"/>\n',strRegion(el),strNum(intensity));
                    fprintf(texfile,'\\color[gray]{%s}\\fill\\%s;\n',strNum(intensity/100),strRegion(el));
                else
                    intensity = max(0,100*intersection{level}.intensity(idel)/maxintens);
                    %fprintf('region: %s, with intensity: %.1f\n',sprintf('%d ',el),intensity/100);
                    id = [];
                    strtex = [];
                    for reg=el(1:end-1)
                        id = sprintf('%sX%s',id,strRegion(reg));
                        strtex = sprintf('%s\\clip\\%s;',strtex,strRegion(reg));
                    end
                    fprintf(svgfile,'\t<use xlink:href="#%s" clip-path="url(#%s)" fill="rgb(%s%%,%s%%,%s%%)"/>\n',strRegion(el(end)),id,strNum(intensity),strNum(intensity),strNum(intensity));
                    %fprintf(svgfile,'\t<use xlink:href="#%s" clip-path="url(#%s)" fill="hsl(0,0%%,%s%%)"/>\n',strRegion(el(end)),id,strNum(intensity)); % does not work with Safari 5.0.4
                    fprintf(texfile,'\\begin{scope}%s\\color[gray]{%s}\\fill\\%s;\\end{scope}\n',strtex,strNum(intensity/100),strRegion(el(end)));
                end
            end
        end
    end
end
%% Footer
fprintf(svgfile,'</g>\n</svg>');
fclose(svgfile);
fprintf(texfile,'\\end{scope}\n\\end{tikzpicture}\n\\end{document}');
fclose(texfile);
fclose(mfile);

%% Compressed file (SVGZ format)
gzip([filename '.svg']);
movefile([filename '.svg.gz'],[filename '.svgz']);

%% Subfunctions
    function str = control2svg(control)
        nodes = control2node(control);
        str = sprintf('M%sQ%s%sT%s',vertex2svg(nodes(1,:)),vertex2svg(control(end,:)),vertex2svg(nodes(2,:)),vertex2svg(nodes([3:end,1],:)));
    end
    function str = control2tex(control)
        nodes = control2node(control);
        str = [];
        Nc = size(control,1);
        for i=1:Nc
            str = sprintf('%s%s..controls%sand%s..',str,vertex2tex(nodes(i,:)),vertex2tex(nodes(i,:)/3+2*control(mod(i-2,Nc)+1,:)/3),vertex2tex(nodes(mod(i,Nc)+1,:)/3+2*control(mod(i-2,Nc)+1,:)/3));
            %str = sprintf('%s\\qbz{%s}{%s}{%s}',str,vertex2tex2(nodes(i,:)),vertex2tex2(control(mod(i-2,Nc)+1,:)),vertex2tex2(nodes(mod(i,Nc)+1,:)));
        end
        str = sprintf('%s%s',str,vertex2tex(nodes(1,:),''));
    end
    function str = vertex2svg(vertex)
        str = [];
        for i=1:size(vertex,1)
            str = sprintf('%s%s,%s ',str,strNum(vertex(i,2)),strNum(vertex(i,1)));
        end
        str= str(1:end);
    end
    function str = vertex2tex(vertex,delimiter)
        if nargin<2, delimiter = '';end
        str = [];
        for i=1:size(vertex,1)
            str = sprintf('%s(%s,%s)%s',str,strNum(vertex(i,2)),strNum(-vertex(i,1)),delimiter);
        end
    end
    function str = vertex2tex2(vertex,delimiter)
        if nargin<2, delimiter = '';end
        str = [];
        for i=1:size(vertex,1)
            str = sprintf('%s%s,%s%s',str,strNum(vertex(i,2)),strNum(-vertex(i,1)),delimiter);
        end
    end
    function str = vertex2m(vertex)
        str = [];
        for i=1:size(vertex,1)
            str = sprintf('%s%s %s;',str,strNum(vertex(i,1)),strNum(vertex(i,2)));
        end
        str = ['[' str(1:end-1) ']'];
    end
    function str = strRegion(idRegion)
        str = sprintf('R');
        idRegion = idRegion-1;
        for d=1:indexDigits
            t = floor(idRegion/10^(indexDigits-d));
            str = sprintf('%s%s',str,char(97+t));
            idRegion = idRegion-10^(indexDigits-d)*t;
        end
    end
    function str = strNum(number)
        if number==0
            str = '0';
        else
            for d = 0:Decimals;
                t = number*10^d;
                if abs(round(t)-t)<10^(d+2)*eps
                    break;
                end
            end
            number = round(t)/10^d;
            if d==0
                str = sprintf('%d',abs(number));
            else
                str = sprintf('%.*f',d,abs(number));
            end
            if abs(number)<1
                str = str(2:end);
            end
            if number<0
                str = ['-' str];
            end
        end
    end
end
