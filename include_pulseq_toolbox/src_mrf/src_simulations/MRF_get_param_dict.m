function P = MRF_get_param_dict(P, constraints)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% ----- inputs: -----
% P: contains the physical parameter ranges and step configs
% constraints: contains a list of physical boundaries

% ----- output: -----
% P: structured dictionary parameters

% -----example: -----
% P.T1.range  = [0.005, 4];   P.T1.factor  = 1.2;
% P.T2.range  = [0.005, 2];   P.T2.factor  = 1.2;
% P.T1p.range = [0.005, 2];   P.T1p.factor = 1.2;
% P.B1.range  = [0.8, 1.2];   P.B1.step    = 0.1;
% P.B0.range  = [-50, 50];    P.B0.step    = 25;
% P.ADC.range = [0, 1.5e-3];  P.ADC.step   = 0.5e-3;
% constraints = {'T2<T1'; 'T2<T1p'; 'T1p<T1'};

p_names    = fieldnames(P);
n_params   = numel(p_names);

% calculate parameter grid
for j = 1:n_params
    temp_name = p_names{j};
    eval(['temp_p = P.' temp_name ';']);
    if isfield(temp_p, 'step')
        temp_p = (temp_p.range(1) : temp_p.step : temp_p.range(2))';
    end
    if isfield(temp_p, 'factor')
        temp_n = floor(log(temp_p.range(2)/temp_p.range(1)) / log(temp_p.factor));
        temp_p = temp_p.range(1) * temp_p.factor.^(0:temp_n)';
    end
    eval([temp_name '= temp_p;']);
end
clear P temp_name temp_p temp_n;

% calculate mesh grids
temp_eval = ['[' p_names{1}];
for j = 2:n_params
    temp_eval = [temp_eval ', ' p_names{j}];
end
temp_eval = [temp_eval '] = ndgrid(' p_names{1}];
for j = 2:n_params
    temp_eval = [temp_eval ', ' p_names{j}];
end
temp_eval = [temp_eval ');'];
eval(temp_eval);
clear temp_eval;

% convert to 1d arrays
for j = 1:n_params
    eval([p_names{j} ' = ' p_names{j} '(:);']);
end

% apply constraints
if nargin>1
    for c = 1:numel(constraints)
        eval(['idx_constraint = ' constraints{c} ';']);
        for j = 1:n_params
            eval([p_names{j} '(~idx_constraint) = [];']);
        end
        clear idx_constraint;
    end
end

% generate output struct
for j = 1:n_params
    eval(['P.' p_names{j} '=' p_names{j} ';']);
end

end

