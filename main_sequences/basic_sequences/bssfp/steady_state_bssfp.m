%% bSSFP Bloch (single isochromat, on-resonance) + analytic steady-state comparison
% Inputs (edit these):
T1    = 1.3;                % [s]
T2    = 0.04;               % [s]
TR    = 9.0e-3;             % [s]
alpha = 60 * pi/180;        % [rad]
NR    = 2000;               % number of readouts
M0    = 1;                  % equilibrium magnetization (proton density scaling)

TE = TR/2;

% Relaxation factors
E1  = exp(-TR/T1);
E2  = exp(-TR/T2);
E1h = exp(-TE/T1);
E2h = exp(-TE/T2);

% State vector [Mx; My; Mz]
M = [0; 0; M0];

Mxy = complex(zeros(NR,1));

%% --- Bloch simulation (single isochromat, on-resonance) ---
for n = 1:NR

    % RF phase cycling for bSSFP: 0, pi, 0, pi, ...
    phi = mod(n-1,2)*pi;

    % Rotation matrices (inline)
    ca = cos(alpha); sa = sin(alpha);
    Rx = [1  0   0;
          0 ca -sa;
          0 sa  ca];

    cp = cos(phi); sp = sin(phi);
    Rz_p  = [ cp -sp 0;
              sp  cp 0;
              0   0  1];
    Rz_m  = [ cp  sp 0;
             -sp  cp 0;
              0   0  1];

    % RF about axis with phase phi: Rz(phi)*Rx(alpha)*Rz(-phi)
    M = Rz_p * Rx * Rz_m * M;

    % Free precession (on-resonance => no z-rotation) + relaxation for TE
    M = [E2h*M(1);
         E2h*M(2);
         E1h*M(3) + (1-E1h)*M0];

    % Readout at TE
    Mxy(n) = M(1) + 1i*M(2);

    % Second half TR: again free precession + relaxation for TE
    M = [E2h*M(1);
         E2h*M(2);
         E1h*M(3) + (1-E1h)*M0];
end

%% --- Analytic steady-state echo amplitude (TE = TR/2, on-resonance) ---
% Scheffler 2003 (Principles and applications of balanced SSFP techniques), Eq. (1): :contentReference[oaicite:1]{index=1}
%   Mss = M0 * sqrt(E2) * (1 - E1) * sin(alpha) / (1 - (E1 - E2)*cos(alpha) - E1*E2)
Mss = M0 * sqrt(E2) * (1 - E1) * sin(alpha) / (1 - (E1 - E2)*cos(alpha) - E1*E2);

%% --- Visualization ---
n = (1:NR).';

figure(1); clf;
plot(n, abs(Mxy), 'LineWidth', 1.5); hold on;
yline(abs(Mss), '--', 'LineWidth', 2.0);
grid on;
xlabel('Readout number');
ylabel('|M_{xy}(TE)|');
title(sprintf('bSSFP single-iso Bloch vs analytic steady state (T1=%.0f ms, T2=%.0f ms, TR=%.1f ms, \\alpha=%.0f°)', ...
    1e3*T1, 1e3*T2, 1e3*TR, alpha*180/pi));
legend('Bloch transient |M_{xy}|', 'Analytic |M_{ss}|', 'Location','best');

figure(2); clf;
plot(n, angle(Mxy), '.');
grid on;
xlabel('Readout number');
ylabel('phase(M_{xy}(TE)) [rad]');
title('bSSFP Bloch phase evolution at TE (with 0/\pi RF phase cycling)');