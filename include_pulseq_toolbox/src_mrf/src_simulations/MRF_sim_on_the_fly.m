function out = MRF_sim_on_the_fly(seq, external_path, flag_mrf)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

out = [];
if flag_mrf == 0
    return;
end

%% write sequence object to .seq file
if exist(external_path, 'file')==0
    [~, external_path] = fileparts(external_path);
    external_path      = ['temp_mrf_' external_path(1:11) '.seq'];
    seq.write(external_path);
end

%% read .seq file
[~, SIM] = MRF_read_seq_file(external_path, 128*1e6, zeros(1000,1), zeros(1000,1), [], 'spiral_out', [], 1e-6, 0);
 
%% dictionary parameters
P.T1  = [5000; 1000; 150] *1e-3;
P.T2  = [ 600;   70;  20] *1e-3;
P.T1p = [ 700;   80;  30] *1e-3;
P.T2p = [ 800;   90;  40] *1e-3;
P.ADC = [0.01; 0.005; 0];
P.T1p_adia = [1200; 600; 300] *1e-3;

%% pre-simulation
SIM = MRF_sim_pre(SIM, P, [], 'EPG', 1, 1);

%% EPG simulation
Mxy = MRF_sim_EPG(SIM, P);

%% EPG simulation for different relax times
figure();
hold on
for j=1:numel(P.T1)
    plot(abs(Mxy(:,j)), '.', 'MarkerSize', 15)   
end
xlabel('number of ADCs')
ylabel('abs(Mxy) [a.u.]')
set(gca, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'Bold', 'LineWidth', 2);
xlim([0 size(Mxy,1)]);
title('Fingerprint: Mxy at start of ADC');
pbaspect([2 1 1])

end

