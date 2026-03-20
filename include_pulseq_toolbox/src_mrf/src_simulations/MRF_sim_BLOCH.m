function M_dict = MRF_sim_BLOCH(SIM, P, z, dw0_distr, comp_energy)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

% this function simulates a dictionary M_dict for different parameters P
% the simulation is based on spin isochromats in z-direction
% rf pulses are simulated as instantaneous events
% slice profiles and relaxation losses during adiabatic pulses are corrected
% the simulation models:
%   T1/T2 relaxation
%   T2' relaxation
%   T1p/T2p relaxation
%   B0 and B1+ deviations
%   apparent diffusion
%   slice profiles
% approximations:
%   pre-simulated slice-profiles for 1:1:270°
%   adiabatic losses approximated via look-up-tables

% ----- input -----
% SIM
%   .ID :  (N_steps x 1) ID for Bloch operator [1...10]
%   .RF :  (N_steps x 1) complex RF waveform [rad/s]
%   .GZ :  (N_steps x 1) z gradient moments [rad/m] or waveforms [rad/s/m]
%   .DB :  (N_steps x 1) b-value of the diffusion encoding tensor [s/mm^2]
%   .DT :  (N_steps x 1) dynamic time steps [s]
%   .PHI:  (N_adc x 1)   phase of receiver [rad]
%   .SPROF: (Niso x 270)  pre-computed slice profiles 1° ... 270°
% P
%   .T1:  (N_dict x 1) T1 relaxation time [s]
%   .T2:  (N_dict x 1) T2 relaxation time [s]
%   .T1p: (N_dict x 1) T1p relaxation time [s]
%   .T2p: (N_dict x 1) T2p relaxation time [s]
%   .m1p: (N_dict x 1) linear T1p dispersion [ms/kHz]
%   .T1p_adia: (N_dict x 1) adiabatic T1p relaxation time [s]
%   .adc  (N_dict x 1) apparent diffusion coefficient (mm^2/s)
%   .db1: (N_dict x 1) B1+ field deviation
%   .dw0: (N_dict x 1) B0 field deviation [rad/s]
%   .TM:  (N x 4x4)    pre-computed transition matrix for adiabatic pulses
% z:           (N_iso x 1) z position of isochromats [m]
% dw0_distr:   (N_iso x 1) B0 field distribution [rad/s]
% comp_energy: 0 -> full dict; >0 -> compressed dict (e.g. 99.99%)

% ----- output -----
% M_dict : (N_adc x N_dict) complex signal

% load pre-compressed dictionary
if nargin<5
    comp_energy = 0; % no compression
end
if comp_energy>0
    [dict_comp, dict_hash] = MRF_find_comp_dict(comp_energy, SIM, P, z, dw0_distr);
    if ~isempty(dict_comp)
        M_dict = dict_comp;
        fprintf(' \n');
        fprintf(['loading compressed dictionary from: ' dict_hash ' \n']);
        fprintf(' \n');
        return;
    end
end

% default parameters
N_dict = numel(P.T1);
if isfield(P, 'T1')
    T1 = P.T1;
else
    error('define T1!');
end
if isfield(P, 'T2')
    T2 = P.T2;
else
    error('define T2!');
end
if isfield(P, 'db1')
    db1 = P.db1;
else
    db1 = ones(N_dict, 1);
end
if isfield(P, 'dw0')
    dw0 = P.dw0;
else
    dw0 = zeros(N_dict, 1);
end
if isfield(P, 'T1p')
    T1p = P.T1p;
else
    T1p = T2;
end
if isfield(P, 'T2p')
    T2p = P.T2p;
else
    T2p = T2;
end
if isfield(P, 'm1p')
    m1p = P.m1p;
else
    m1p = zeros(N_dict, 1);
end
if isfield(P, 'T1p_adia')
    T1p_adia = P.T1p_adia;
else
    T1p_adia = T1p;
end
if isfield(P, 'ADC')
    ADC = P.ADC;
else
    ADC = zeros(N_dict, 1);
end
clear P;

% operator lists
ID    = SIM.ID;
RF    = SIM.RF;
GZ    = SIM.GZ;
DB    = SIM.DB;
DT    = SIM.DT;
PHI   = SIM.PHI;
TM    = SIM.TM;
SPROF = SIM.SPROF;
clear SIM;

% init
N_steps = numel(ID);
N_iso   = numel(z);
N_adc   = sum(ID==0);
M_dict  = zeros(N_adc, N_dict);

% add dummy dimensions to z
z = reshape(z, 1, 1, []);

tic;
fprintf('\n');
fprintf(['calculating dictionary: ' num2str(N_dict) ' atoms via Bloch \n']);

% calculate dictionary
parfor j = 1:N_dict

    % init magnetization, Mz = 1
    M = zeros(4, N_iso);
    M(3:4,:) = 1;
    M_sim = zeros(N_adc, 1);

    % B1+ deviation
    RF_ = RF .* (db1(j)*(ID~=3) + (ID==3)); % apply for all RF events, except adiabatic pulses

    % B0 deviation
    if ~isempty(dw0_distr)
        dw0_ = reshape(dw0(j) + dw0_distr, 1, 1, []);
    else
        dw0_ = dw0(j);
    end

    % counter for ADCs
    n_adc = 1;
    
    % forward simulation using different Bloch operators
    for k = 1:N_steps
    switch(ID(k))
    
        case 0 % adc
            M_sim(n_adc,1) = mean(M(1,:) + 1i * M(2,:)) * exp(-1i * PHI(n_adc));
            n_adc = n_adc + 1;                

        case 1 % free relaxation & dephasing
            M = sim_free_relax_deph(M, DT(k), T1(j), T2(j), DT(k)*dw0_);

        case 2 % global rf excitation
            M = sim_rf(M, abs(RF_(k)), angle(RF_(k)));

        case 3 % adiabatic rf excitation
            M = squeeze(TM(j,abs(RF_(k)),:,:)) * M;

        case 4 % spin-lock
            M = sim_spin_lock(M, DT(k), RF_(k), T1p(j), T2p(j), m1p(j), dw0_);

        case {5, 6} % slice selective rf excitation or refocusing
            alpha = abs(RF_(k));
            phi   = angle(RF_(k));
            id    = max([1, round(alpha*180/pi)]);
            alpha = alpha * abs(SPROF(:,id));
            phi   = phi   + angle(SPROF(:,id));
            M     = sim_rf(M, alpha, phi);

        case {8, 9} % perfect crusher
            M(1,:) = 0;
            M(2,:) = 0;

        case 10 % z gradient
            M = sim_z_deph(M, GZ(k)*z);

        case 11 % diffusion
            M = sim_diffusion(M, DB(k), ADC(j));

        case 12 % adiabatic spin-lock
            M = sim_rf(M, pi, 0);
            M = exp(-DT(k)/T1p_adia(j)) * M;
            M(4,:) = 1;

    end
    end

    M_dict(:,j) = M_sim;

end

temp_t = toc;
fprintf(['   ' num2str(temp_t, '%.1f') 's ... complete! '  num2str(temp_t/N_dict*1e3, '%.1f')  ' ms/atom \n']);

% save pre-compressed dictionary
if comp_energy>0
    [~, S ,V] = svd(M_dict.', 'econ');
    svals     = diag(S)/S(1);
    svals     = cumsum(svals.^2./sum(svals.^2));
    NPCs      = find(svals > comp_energy, 1, 'first');
    dict_phi  = single(V(:,1:NPCs));
    dict_svd  = (M_dict.'*dict_phi).';
    dict_comp.phi  = dict_phi;
    dict_comp.svd  = dict_svd;
    dict_comp.NPCs = NPCs;
    save(dict_hash, 'dict_comp', '-v7.3');
    M_dict = dict_comp;
end

end

%% ------------------- simulate: depashing -------------------
function M = sim_z_deph(M, psi)
    cos_psi = cos(psi);
    sin_psi = sin(psi);
    O       = zeros(size(psi));
    I       = ones(size(psi));    
    B       = [ cos_psi, sin_psi, O, O;
               -sin_psi, cos_psi, O, O;
                O,       O,       I, O;
                O,       O,       O, I ];
    M = matrix_vector_iso_prod(M, B);
end

%% ---------- simulate: free relaxation & depashing ----------
function M = sim_free_relax_deph(M, dt, T1, T2, psi)
    E1      = exp(-dt/T1);
    E2      = exp(-dt/T2);
    cos_psi = cos(psi);
    sin_psi = sin(psi);
    O       = zeros(size(psi));
    I       = ones(size(psi));
    B       = [ E2*cos_psi, E2*sin_psi, O,    O;
               -E2*sin_psi, E2*cos_psi, O,    O;
                O,          O,          I*E1, I-E1;
                O,          O,          O,    I ];
    M = matrix_vector_iso_prod(M, B);
end

%% ---------- simulate: rf -----------
function M = sim_rf(M, alpha, phi)
    if numel(alpha)>1
        alpha = reshape(alpha,   1, 1, []);
        phi   = reshape(phi,   1, 1, []);
    end
    cos_phi   = cos(phi);
    sin_phi   = sin(phi);
    cos_alpha = cos(alpha);
    sin_alpha = sin(alpha);
    O         = zeros(size(alpha));
    I         = ones(size(alpha));
    B = [ cos_phi.^2 + sin_phi.^2 .* cos_alpha,   cos_phi .* sin_phi .* (1- cos_alpha),   -sin_phi .* sin_alpha,  O;
          cos_phi .* sin_phi .* (1- cos_alpha),   sin_phi.^2 + cos_phi.^2 .* cos_alpha,    cos_phi .* sin_alpha,  O;
          sin_phi .* sin_alpha,                  -cos_phi .* sin_alpha,                    cos_alpha,             O;
          O,                                      O,                                       O,                     I ];
    M = matrix_vector_iso_prod(M, B);
end

%% ---------- simulate: diffusion ----------
function M = sim_diffusion(M, B, ADC)
    B = [ exp(-B*ADC), 0, 0, 0;
          0, exp(-B*ADC), 0, 0;
          0, 0,           1, 0;
          0, 0,           0, 1 ];
    M = B * M;
end

%% ------------------- simulate: spin-lock -------------------
function M = sim_spin_lock(M, tSL, wSL, T1p, T2p, m1p, dw0)

    % assume linear T1p dispersion
    T1p = T1p + m1p * abs(wSL) * 1e-6 / 2 / pi; % m1p -> ms/kHz
    
    % dummy dimensions for isochromats
    I = ones(size(dw0));

    % rotating frame relaxation depending on SL phase
    phi = angle(wSL);
    R1p = 1/T1p;
    R2p = 1/T2p;
    Rx  = (R1p * cos(phi)^2 + R2p * sin(phi)^2) * I;
    Ry  = (R1p * sin(phi)^2 + R2p * cos(phi)^2) * I;
    Rz  = R2p * I;
    wSL = wSL * I;
 
    % define bloch operator for spin-locking
    BSL = [ -Rx,         dw0,        -imag(wSL);
            -dw0,       -Ry,          real(wSL);
            imag(wSL),  -real(wSL),  -Rz ];
    
    % calculate matrix exponential and switch to 4x4
    B = zeros(3, 3, numel(dw0));    
    for j=1:numel(dw0)
        B(:,:,j) = expm(squeeze(BSL(:,:,j))*tSL);
    end
    B(4,4,:) = 1;
    M = matrix_vector_iso_prod(M, B);

end

%% ------------------- sliced matrix vector product -------------------
function M_new = matrix_vector_iso_prod(M, B)

    % this function is used for "sliced" matrix vector products
    % the Bloch matrix B can be different for each individual isochromat (4x4xn)
    % in this case, the matrix vector multiplication needs to be "sliced"
    % if B is a 4x4 matrix, the function is reduced to a simple matrix vector product

    % ----- input: -----
    % M: 4 x n
    % B: 4 x 4 or 4 x 4 x n
    
    % ----- output: -----
    % M_new: 4 x n  

    if ismatrix(B)
        M_new = B * M;
    else

        % C++ mex version
        M_new = mex_matrix_vector_iso_prod(M, B);

        % Matlab version:
        % n = size(M,2);
        % M_new = zeros(4, n);
        % for j=1:n
        %     M_new(1,j) = B(1,1,j)*M(1,j) + B(1,2,j)*M(2,j) + B(1,3,j)*M(3,j) + B(1,4,j)*M(4,j);
        %     M_new(2,j) = B(2,1,j)*M(1,j) + B(2,2,j)*M(2,j) + B(2,3,j)*M(3,j) + B(2,4,j)*M(4,j);
        %     M_new(3,j) = B(3,1,j)*M(1,j) + B(3,2,j)*M(2,j) + B(3,3,j)*M(3,j) + B(3,4,j)*M(4,j);
        %     M_new(4,j) = B(4,1,j)*M(1,j) + B(4,2,j)*M(2,j) + B(4,3,j)*M(3,j) + B(4,4,j)*M(4,j);
        % end

    end

end

% ------------------- mex_matrix_vector_iso_prod.cpp -------------------
% #include "mex.h"
% 
% void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
% {
%     // input:
%     const double *M = mxGetPr(prhs[0]); // M: 4 x n
%     const double *B = mxGetPr(prhs[1]); // B: 4 x 4 x n
%     mwSize n = mxGetN(prhs[0]);         // number of isochromats
% 
%     // output: 4 x n
%     mwSize dims[2] = {4, n};
%     plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
%     double *M_new = mxGetPr(plhs[0]);
% 
%     // sliced matrix vector product
%     for (mwSize j = 0; j < n; ++j)
%     {
%         const double *M_j = M     + 4  * j;  // M(:,j)
%         const double *B_j = B     + 16 * j;  // B(:,:,j)
%         double *M_new_j   = M_new + 4  * j;  // M_new(:,j)
% 
%         M_new_j[0] = B_j[0] * M_j[0] + B_j[4] * M_j[1] + B_j[8]  * M_j[2] + B_j[12] * M_j[3];
%         M_new_j[1] = B_j[1] * M_j[0] + B_j[5] * M_j[1] + B_j[9]  * M_j[2] + B_j[13] * M_j[3];
%         M_new_j[2] = B_j[2] * M_j[0] + B_j[6] * M_j[1] + B_j[10] * M_j[2] + B_j[14] * M_j[3];
%         M_new_j[3] = B_j[3] * M_j[0] + B_j[7] * M_j[1] + B_j[11] * M_j[2] + B_j[15] * M_j[3];
%     }
% }