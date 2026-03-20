% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

%% try mex compilation
eval(['mex -I' [pulseq_get_path('find_eigen_path') 'eigen-3.4.0/'] ' mex_BLOCH_euler.cpp']);
eval(['mex -I' [pulseq_get_path('find_eigen_path') 'eigen-3.4.0/'] ' mex_BLOCH_rk4.cpp']);
eval(['mex -I' [pulseq_get_path('find_eigen_path') 'eigen-3.4.0/'] ' mex_BLOCH_expm.cpp']);
eval(['mex -I' [pulseq_get_path('find_eigen_path') 'eigen-3.4.0/'] ' mex_BLOCH_SL_euler.cpp']);
eval(['mex -I' [pulseq_get_path('find_eigen_path') 'eigen-3.4.0/'] ' mex_BLOCH_SL_rk4.cpp']);
eval(['mex -I' [pulseq_get_path('find_eigen_path') 'eigen-3.4.0/'] ' mex_BLOCH_SL_expm.cpp']);
eval(['mex -I' [pulseq_get_path('find_eigen_path') 'eigen-3.4.0/'] ' mex_EPG_rf_relax.cpp']);
mex mex_matrix_vector_iso_prod.cpp;

% compiled on Win 11: Building with 'Microsoft Visual C++ 2022'
% compiled on Ubuntu: Building with 'g++'