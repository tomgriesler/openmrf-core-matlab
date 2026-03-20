%% DefineRect.m
%
% SEE: phantom
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 19-05-2011 (dd-mm-yyyy)

angle = pi/3;
R = [cos(angle) -sin(angle);sin(angle) cos(angle)];
width = 0.5*[1,0.8];
rect.FOV = [1,1];
rect.region = cell(1,1);
rect.region{1} = struct( 'type', 'polygon', 'weight', [1], 'vertex', [1 -1; 1 1; -1 1;-1 -1]*diag(width/2)*R);
%clear angle R width;