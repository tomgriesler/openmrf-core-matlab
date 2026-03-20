// -------------------------------------------------------------
// Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
// -------------------------------------------------------------

// ---------------------- use in Matlab: -----------------------
// T = mex_EPG_rf_relax(w1x, w1y, R1, R2, dt)
// w1x:  N x 1 double array
// w1y:  N x 1 double array
// R1:   1 x 1 double
// R2:   1 x 1 double
// dt:   1 x 1 double
// -------------------------------------------------------------

// -------------------------------------------------------------
// -------------------------- MEX, C++ -------------------------
// -------------------------------------------------------------

#include <mex.h>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
# define PI 3.14159265358979323846264338
using namespace std;
using namespace Eigen;

// dimension of Bloch matrix: 3x3
constexpr int n = 3;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{

    // define pointers to the input arrays
    double* w1x = mxGetPr(prhs[0]);     // Nx1 double array
    double* w1y = mxGetPr(prhs[1]);     // Nx1 double array
    double R1   = mxGetScalar(prhs[2]); // 1x1 double scalar
    double R2   = mxGetScalar(prhs[3]); // 1x1 double scalar
    double dt   = mxGetScalar(prhs[4]); // 1x1 double scalar

    // get number of dt steps
    int N = mxGetNumberOfElements(prhs[0]);

    // prepare Eigen arrays
    MatrixXcd B_RF(3, 3);
    MatrixXcd B_Relax(3, 3);
    MatrixXcd T(3, 3);
    B_RF.setZero();
    B_Relax.setZero();
    T.setZero();
    T(0,0) = std::complex<double>(1.0, 0.0);
    T(1,1) = std::complex<double>(1.0, 0.0);
    T(2,2) = std::complex<double>(1.0, 0.0);   

    // Compute Transition Matrix
    for (int j = 0; j < N; ++j)
    {
        double alpha = sqrt(w1x[j]*w1x[j] + w1y[j]*w1y[j]) * dt;
        double phi   = atan2(w1y[j], w1x[j]);
        B_Relax << exp(-dt*R2),    0.0,            0.0,
                   0.0,            exp(-dt*R2),    0.0,
                   0.0,            0.0,            exp(-dt*R1);
        B_RF << cos(alpha/2.0)*cos(alpha/2.0),                                          exp(std::complex<double>(0, 2.0*phi))*sin(alpha/2.0)*sin(alpha/2.0),    exp(std::complex<double>(0, phi-PI/2.0))*sin(alpha),
                exp(std::complex<double>(0, -2.0*phi))*sin(alpha/2.0)*sin(alpha/2.0),   cos(alpha/2.0)*cos(alpha/2.0),                                          exp(std::complex<double>(0, PI/2.0-phi))*sin(alpha),
                exp(std::complex<double>(0, -phi-PI/2.0))*sin(alpha)/2.0,               exp(std::complex<double>(0, phi+PI/2.0))*sin(alpha)/2.0,                cos(alpha);
        
        T = B_RF * B_Relax * T;
    }

    // output final transition matrix
    plhs[0] = mxCreateDoubleMatrix(n, n, mxCOMPLEX);
    double *realT = mxGetPr(plhs[0]);
    double *imagT = mxGetPi(plhs[0]);
    for (int j = 0; j < 9; ++j)
    {
        realT[j] = std::real(T(j));
        imagT[j] = std::imag(T(j));
    }

}