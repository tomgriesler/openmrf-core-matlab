// -------------------------------------------------------------
// Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
// -------------------------------------------------------------

// ---------------------- use in Matlab: -----------------------
// T = mex_BLOCH_SL_expm(M, w1x, w1y, dw0, R1p, R2p, dt)
// M:    4 x 1 double array
// w1x:  N x 1 double array
// w1y:  N x 1 double array
// dw0:  N x 1 double array
// R1p:  1 x 1 double
// R2p:  1 x 1 double
// dt:   1 x 1 double
// -------------------------------------------------------------

// -------------------------------------------------------------
// -------------------------- MEX, C++ -------------------------
// -------------------------------------------------------------

#include <mex.h>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
using namespace std;
using namespace Eigen;

// dimension of Bloch matrix: 3x3
constexpr int n = 3;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{

    // define pointers to the input arrays
    double* Minit = mxGetPr(prhs[0]);     // 4x1 double array
    double* w1x   = mxGetPr(prhs[1]);     // Nx1 double array
    double* w1y   = mxGetPr(prhs[2]);     // Nx1 double array
    double* dw0   = mxGetPr(prhs[3]);     // Nx1 double array
    double R1p    = mxGetScalar(prhs[4]); // 1x1 double scalar
    double R2p    = mxGetScalar(prhs[5]); // 1x1 double scalar
    double dt     = mxGetScalar(prhs[6]); // 1x1 double scalar

    // get number of dt steps
    int N = mxGetNumberOfElements(prhs[1]);

    // prepare Eigen arrays
    Eigen::Matrix<double, n, 1> M;
    Eigen::Matrix<double, n, n> B;
    Eigen::Matrix<double, n, n> T;
    T.setIdentity();

    // copy initial M into Eigen vector
    for (int j = 0; j < n; ++j)
    M(j) = Minit[j];

    // Compute Transition Matrix
    for (int j = 0; j < N; ++j)
    {
        B << -R1p,      dw0[j],  -w1y[j],
             -dw0[j],  -R2p,      w1x[j],
              w1y[j],  -w1x[j],  -R2p;
        B = B * dt;
        T = B.exp() * T;
    }

    // Compute the final magnetization
    M = T * M;

    // output final M
    plhs[0] = mxCreateDoubleMatrix(4, 1, mxREAL);
    double* output_M = mxGetPr(plhs[0]);
    for (int j = 0; j < n; ++j)
    output_M[j] = M(j);
    output_M[3] = 1.0;

    // output final transition matrix
    plhs[1] = mxCreateDoubleMatrix(4, 4, mxREAL);
    double* output_T = mxGetPr(plhs[1]);
    output_T[0]  = T(0,0);
    output_T[1]  = T(1,0);
    output_T[2]  = T(2,0);
    output_T[3]  = 0.0;
    output_T[4]  = T(0,1);
    output_T[5]  = T(1,1);
    output_T[6]  = T(2,1);
    output_T[7]  = 0.0;
    output_T[8]  = T(0,2);
    output_T[9]  = T(1,2);
    output_T[10] = T(2,2);
    output_T[11] = 0.0;
    output_T[12] = 0.0;
    output_T[13] = 0.0;
    output_T[14] = 0.0;
    output_T[15] = 1.0;

}