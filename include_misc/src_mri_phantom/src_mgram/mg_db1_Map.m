function B1_Map = mg_db1_Map(Nxy, B1max, B1min, mask, radius)

if isempty(mask)
    mask = ones(Nxy,Nxy);
end
if isempty(radius)
    radius = 0.95;
end

B1_Map = mg_get_tukey_window(Nxy, Nxy, radius, 0)  .* mask;
B1_Map = (B1_Map - min(B1_Map(mask==1))) .* mask;
B1_Map = B1_Map / max(B1_Map(mask==1));
B1_Map = (B1_Map * (B1max-B1min) + B1min) .* mask;

end

