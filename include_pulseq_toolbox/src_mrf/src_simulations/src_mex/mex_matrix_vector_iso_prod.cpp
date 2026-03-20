// -------------------------------------------------------------
// Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
// -------------------------------------------------------------

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // input:
    const double *M = mxGetPr(prhs[0]); // M: 4 x n
    const double *B = mxGetPr(prhs[1]); // B: 4 x 4 x n
    mwSize n = mxGetN(prhs[0]);         // number of isochromats

    // output: 4 x n
    mwSize dims[2] = {4, n};
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    double *M_new = mxGetPr(plhs[0]);

    // sliced matrix vector product
    for (mwSize j = 0; j < n; ++j)
    {
        const double *M_j = M     + 4  * j;  // M(:,j)
        const double *B_j = B     + 16 * j;  // B(:,:,j)
        double *M_new_j   = M_new + 4  * j;  // M_new(:,j)

        M_new_j[0] = B_j[0] * M_j[0] + B_j[4] * M_j[1] + B_j[8]  * M_j[2] + B_j[12] * M_j[3];
        M_new_j[1] = B_j[1] * M_j[0] + B_j[5] * M_j[1] + B_j[9]  * M_j[2] + B_j[13] * M_j[3];
        M_new_j[2] = B_j[2] * M_j[0] + B_j[6] * M_j[1] + B_j[10] * M_j[2] + B_j[14] * M_j[3];
        M_new_j[3] = B_j[3] * M_j[0] + B_j[7] * M_j[1] + B_j[11] * M_j[2] + B_j[15] * M_j[3];
    }
}