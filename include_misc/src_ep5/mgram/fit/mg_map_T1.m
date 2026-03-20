function [T1_Map, M0_Map, Eff_Map, R2_Map, mask_fit] = mg_map_T1(ImageArr, TI, mask_fit)

% V1: Maximilian Gram; 21.03.2024;
% V2: Maximilian Gram; 06.06.2024;
%     fit model accoring to doi.org/10.1002/mrm.28779
%     S = P(1) * (1-(1+P(2))*exp(-TI/T1))

if nargin<3
    mask_fit = mg_get_mask_fit(squeeze(mean(ImageArr)));
end
if size(TI,1) < size(TI,2)
    TI = TI';
end

%% pixel:pixel fit

Ny = size(ImageArr,2);
Nx = size(ImageArr,3);

T1_Map    = zeros(Ny,Nx);
M0_Map    = zeros(Ny,Nx);
Eff_Map   = zeros(Ny,Nx);
R2_Map    = zeros(Ny,Nx);

parfor y = 1:Ny
   for x = 1:Nx
   if mask_fit(y,x)==1

       data    = squeeze(ImageArr(:,y,x));
       P_start = [abs(data(end)), 1.0, 1.0];

       opts   = optimset('Display','off');       
       FITfun = @(P) P(1) * (1-(1+P(2))*exp(-TI/P(3)));
       RSSfun = @(P) sum((data - FITfun(P)).^2);
       P_fit  = fminsearch(RSSfun, P_start, opts);

       RSS    = RSSfun(P_fit);
       TSS    = sum( (data-mean(data)).^2 );
       R2     = 1 - RSS/TSS;

       M0_Map(y,x)  = P_fit(1); 
       Eff_Map(y,x) = P_fit(2);
       T1_Map(y,x)  = P_fit(3);
       R2_Map(y,x)  = R2;
       
   end
   end
end

end

