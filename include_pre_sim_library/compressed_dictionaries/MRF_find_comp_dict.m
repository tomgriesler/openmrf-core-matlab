function [dict_comp, dict_hash] = MRF_find_comp_dict(comp_energy, SIM, P, z, dw0_distr)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% combine all inputs in a single 1d array
inputs_all    = [];
temp_fields = fieldnames(SIM);
for k = 1:numel(temp_fields)
    inputs_all = [inputs_all; single(SIM.(temp_fields{k})(:))];
end
temp_fields = fieldnames(P);
for k = 1:numel(temp_fields)
    inputs_all = [inputs_all; single(P.(temp_fields{k})(:))];
end
if nargin > 3
    inputs_all = [inputs_all; single(z)];
end
if nargin > 4
    inputs_all = [inputs_all; single(dw0_distr)];
end
clear temp_fields;

% create hash
dict_hash = [pulseq_get_path('MRF_find_comp_dict') SIM.seq_name '_' pulseq_get_wave_hash([real(inputs_all); imag(inputs_all)]) strrep(num2str(comp_energy),'0.','_') '.mat'];

% find compressed dictionary   
if isfile(dict_hash)
     load(dict_hash);
else
    dict_comp = [];
end

end