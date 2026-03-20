function phantom = mg_phantom_properties(phantom, varargin)

% add properties for all inputs
for j=1:numel(varargin)
    prop = varargin{j};
    if numel(prop)~=phantom.n_phant+1
        error(['wrong number in: ' inputname(j+1)]);
    end
    prop2D = zeros(phantom.Nxy, phantom.Nxy);
    for p = 1 : phantom.n_phant+1        
        prop2D(phantom.ind_map==p) = prop(p);
    end
    prop1D = prop2D(phantom.mask2D);
    eval(['phantom.p.'   inputname(j+1) '=prop;']);
    eval(['phantom.p2D.' inputname(j+1) '=prop2D;']);
    eval(['phantom.p1D.' inputname(j+1) '=prop1D;']);
    clear prop1D prop2D;
end

% defaults for T2*, df0 and db1
if ~isfield(phantom.p, 'T2s')
    phantom.p.T2s   = phantom.p.PD*0   + 1e12;
    phantom.p1D.T2s = phantom.p1D.PD*0 + 1e12;
    phantom.p2D.T2s = phantom.p2D.PD*0 + 1e12;
end
if ~isfield(phantom.p, 'dw0')
    phantom.p.dw0   = phantom.p.PD*0;
    phantom.p1D.dw0 = phantom.p1D.PD*0;
    phantom.p2D.dw0 = phantom.p2D.PD*0;
end
if ~isfield(phantom.p, 'db1')
    phantom.p.db1   = phantom.p.PD*0   + 1;
    phantom.p1D.db1 = phantom.p1D.PD*0 + 1;
    phantom.p2D.db1 = phantom.p2D.PD*0 + 1;
end

% add property noise
if isfield(phantom, 'noise_level')
    rng("default")
    temp_names = fieldnames(phantom.p);
    for j = 1:numel(temp_names)        
        eval(['phantom.p1D.' temp_names{j} ' = phantom.p1D.' temp_names{j} '.* (1 + randn(size(phantom.p1D.' temp_names{j} ')) * phantom.noise_level / 100);']);
        eval(['phantom.p2D.' temp_names{j} '(phantom.pos_map(phantom.mask2D))  = phantom.p1D.' temp_names{j} ';']);
    end
else
    phantom.noise_level = 0;
end

% add masks and indices
phantom.p1D.ind   = phantom.ind_map(phantom.mask2D);
phantom.p2D.ind   = phantom.ind_map;
phantom.p1D.pos   = phantom.pos_map(phantom.mask2D);
phantom.p2D.pos   = phantom.pos_map;
phantom.p1D.mask  = phantom.mask1D;
phantom.p2D.mask  = phantom.mask2D;
phantom.p1D.cmaps = phantom.cmaps(:,phantom.mask1D).';
phantom.p2D.cmaps = phantom.cmaps;

phantom = rmfield(phantom, {'cmaps'; 'ind_map'; 'pos_map'; 'mask2D'; 'mask1D'});

end

