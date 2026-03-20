// -------------------------------------------------------------
// Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
// -------------------------------------------------------------

// ---------------------- use in Matlab: -----------------------
// M = mex_BLOCH_SL_euler(M, w1x, w1y, dw0, R1p, R2p, dt)
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
    double* Minit = mxGetPr(prhs[0]);     // 3x1 double array
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

    // copy initial M into Eigen vector
    for (int j = 0; j < n; ++j)
    M(j) = Minit[j];

    // Euler Integration
    for (int j = 0; j < N; ++j)
    {
        B << -R1p,      dw0[j],  -w1y[j],
             -dw0[j],  -R2p,      w1x[j],
              w1y[j],  -w1x[j],  -R2p;
        M += dt * B * M;
    }

    // output final M
    plhs[0] = mxCreateDoubleMatrix(4, 1, mxREAL);
    double* output_M = mxGetPr(plhs[0]);
    for (int j = 0; j < n; ++j)
    output_M[j] = M(j);
    output_M[3] = 1.0;

}