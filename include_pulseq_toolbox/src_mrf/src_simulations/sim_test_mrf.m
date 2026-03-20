clear

% seq_file = 'test_mrf.seq';
seq_file = 'test_cmrf.seq';

[SEQ, SIM] = MRF_read_seq_file( seq_file, ... % .seq file path
                                128*1e6, ...  % larmor frequency
                                [], ...       % adc time stamps
                                [], ...       % soft delays
                                [], ...       % kz part
                                [], ...       % echo mode
                                1e-6, ...     % simulation raster time
                                0);           % plot flag

%% params

% slice & isochromats
Niso = 1000;
dz   = 8*1e-3;
z    = linspace(-1/2, 1/2, Niso)' * 2 * dz;

% dictionary parameters
N_dict = 14;
[P.T1, P.T2] = NIST_references('3.0T_MnCl2');
P.T1   = linspace(0.1, 2, N_dict)';
P.T1p  = P.T2 .* linspace(1.0, 1.4, N_dict)';
P.T2p  = P.T2 .* linspace(1.4, 1.8, N_dict)';
P.db1  = 0.1  * linspace(-1, 1, N_dict)' + 1;
P.dw0  = 10   * linspace(-1, 1, N_dict)' * 2*pi;

% T2' effect
dw0_distr = randn(Niso, 1) * 10 * 2*pi;

%% pre-simulations
SIM_BLOCH = MRF_sim_pre(SIM, P, z, 'BLOCH', 1, 1);
SIM_EPG   = MRF_sim_pre(SIM, P, z, 'EPG',   1, 1);

%% brute force Bloch simulation
M_dict1 = MRF_sim_brute_force(SIM, P, z, dw0_distr);

%% fast Bloch simulation
M_dict2 = MRF_sim_BLOCH(SIM_BLOCH, P, z, dw0_distr);

%% fast EPG simulation
M_dict3 = MRF_sim_EPG(SIM_EPG, P);

%% vis: compare all

mylim = max([max(abs(real([M_dict1, M_dict2, M_dict3]))), max(abs(real([M_dict1, M_dict2, M_dict3])))]);
mylim = [-1, 1] * mylim;

cmp_real = winter(N_dict);
cmp_imag = autumn(N_dict);

figure();
ax1 = subplot(3,2,1);
hold on
for j=1:N_dict
    plot(real(M_dict1(:,j)), '-', 'LineWidth', 1, 'Color', cmp_real(j,:))
end
ylim(mylim)
title('real(brute force)')

ax2 = subplot(3,2,2);
hold on
for j=1:N_dict
    plot(imag(M_dict1(:,j)), '-', 'LineWidth', 1, 'Color', cmp_imag(j,:))
end
ylim(mylim)
ylim(mylim)
title('imag(brute force)')

ax3 = subplot(3,2,3);
hold on
for j=1:N_dict
    plot(real(M_dict2(:,j)), '-', 'LineWidth', 1, 'Color', cmp_real(j,:))
end
ylim(mylim)
ylim(mylim)
title('real(bloch)')

ax4 = subplot(3,2,4);
hold on
for j=1:N_dict
    plot(imag(M_dict2(:,j)), '-', 'LineWidth', 1, 'Color', cmp_imag(j,:))
end
ylim(mylim)
ylim(mylim)
title('imag(bloch)')

ax5 = subplot(3,2,5);
hold on
for j=1:N_dict
    plot(real(M_dict3(:,j)), '-', 'LineWidth', 1, 'Color', cmp_real(j,:))
end
ylim(mylim)
ylim(mylim)
title('real(epg)')

ax6 = subplot(3,2,6);
hold on
for j=1:N_dict
    plot(imag(M_dict3(:,j)), '-', 'LineWidth', 1, 'Color', cmp_imag(j,:))
end
ylim(mylim)
ylim(mylim)
title('imag(epg)')

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6], 'xy')

%% vis: compare burte fore and fast bloch

mylim = max([max(abs(real([M_dict1-M_dict2]))), max(abs(real([M_dict1-M_dict2])))]);
mylim = [-1, 1] * mylim;

figure();
ax1 = subplot(1,2,1);
hold on
for j=1:N_dict
    plot(real(M_dict2(:,j)-M_dict1(:,j)), '-', 'LineWidth', 1, 'Color', cmp_real(j,:))
end
ylim(mylim)
title('real(brute force - fast bloch)')

ax2 = subplot(1,2,2);
hold on
for j=1:N_dict
    plot(imag(M_dict2(:,j)-M_dict1(:,j)), '-', 'LineWidth', 1, 'Color', cmp_imag(j,:))
end
ylim(mylim)
ylim(mylim)
title('imag(brute force - fast bloch)')

linkaxes([ax1 ax2], 'xy')

%% vis bloch vs epg: ideal slice profile

P_ideal.T1  = 0.8;
P_ideal.T2  = 0.05;
P_ideal.T1p = 0.08;
P_ideal.T2p = 0.08;
P_ideal.dw0 = 0;
P_ideal.db1 = 1;

SIM_BLOCH_ideal = MRF_sim_pre(SIM, P_ideal, z, 'BLOCH', 0, 0);
SIM_EPG_ideal   = MRF_sim_pre(SIM, P_ideal, z, 'EPG',   0, 0);
SIM_BLOCH_ideal.SPROF = SIM_BLOCH_ideal.SPROF*0 + 1;

M_ideal_bloch = MRF_sim_BLOCH(SIM_BLOCH_ideal, P_ideal, z, []);
M_ideal_epg   = MRF_sim_EPG(SIM_EPG_ideal,     P_ideal);

figure()
subplot(1,2,1)
hold on
plot(real(M_ideal_bloch), 'b-', 'LineWidth', 2)
plot(real(M_ideal_epg),   'r-', 'LineWidth', 2)
title('real()')
subplot(1,2,2)
hold on
plot(imag(M_ideal_bloch), 'b-', 'LineWidth', 2)
plot(imag(M_ideal_epg),   'r-', 'LineWidth', 2)
title('imag()')
legend('bloch', 'epg')
