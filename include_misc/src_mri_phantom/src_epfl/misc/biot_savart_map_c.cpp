// biot_savart_map_c.cpp
//
// Computes the magnetic field generated on a Cartesian grid by a circular coil.
//
// 
// INPUTS:	x0      vector of x-coordinates of the map samples (M its length)
//          y0      vector of y-coordinates of the map samples (N its length)
//          R       radius of the coil
//          D       distance between the center of the sample and the center of the coil
//          theta   coil's angular position in the sample frame of reference
//
// OUTPUT:  Map of magnetic field complex values (real part: x-component; imaginary part: y-component)
//
// This is a MEX-file for MATLAB.  
//
// Matthieu Guerquin-Kern, Biomedical Imaging Group - EPF Lausanne, 2008-02-28
// revised in october 2010

#include <stdlib.h>
#include <math.h>
#include "mex.h"

#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#define	MAX(A, B)	((A) > (B) ? (A) : (B))

/* Input Arguments */

#define	x0_IN       prhs[0]
#define	y0_IN       prhs[1]
#define	R_IN        prhs[2]
#define	D_IN        prhs[3]
#define	theta_IN	prhs[4]

/* Output Arguments */

#define	S_OUT   plhs[0]

// Global variables
static	double	mu = 1e-7;		// vacuum magnetic permitivity divided by 2*pi
static  int     Nb_angles = 40; // number of sides of the polygon that approximates the coil (the higher, the better)
static  double  I  =  1;		// Intensity in the coil (reciprocity principle)

double norm( double X[]) // Computes the norm of a vector of 3 scalars
{
    return sqrt(pow(X[0],2)+pow(X[1],2)+pow(X[2],2));
}

void cross( double X[], double Y[], double Z[]) // computes the cross-product of two vectors of 3 scalars
{
    Z[0] = X[1]*Y[2]-X[2]*Y[1];
    Z[1] = X[2]*Y[0]-X[0]*Y[2];
    Z[2] = X[0]*Y[1]-X[1]*Y[0];
    return;
}

static void biot_savart( double B[], const double M[], const double RC[], const double RS[], const double I )
{
    double	D_theta,cst;
    double S[]={0,0,0},dL[]={0,0,0},C[]={0,0,0};
	int theta;
        
    B[0] = 0;
    B[1] = 0;
    B[2] = 0;
    D_theta = 2.0*M_PI/(double)Nb_angles;
    
    for (theta=0;theta<Nb_angles;theta++){
        S[0] = M[0] - RC[theta];
        S[1] = M[1] - RS[theta];
        S[2] = M[2];
        dL[0] =  D_theta * RS[theta];
        dL[1] = -D_theta * RC[theta];
        dL[2] = 0;
        cross(dL,S,C);
        cst = mu*I/pow(norm(S),3);
        B[0] = B[0] + cst*C[0];
        B[1] = B[1] + cst*C[1];
        B[2] = B[2] + cst*C[2];
    }
    return;
}

void rotation(double M[],const double c,const double s){
    double x=M[0],y=M[1];
    
    M[0] =  c * x + s * y;
    M[1] = -s * x + c * y;
}

void permutation(double M[]){
    double x= M[0],y= M[1],z= M[2];
    
    M[0] = z;
    M[1] = x;
    M[2] = y;
}

void inv_permutation(double M[]){
    double x= M[0],y= M[1],z= M[2];
    
    M[0] = y;
    M[1] = z;
    M[2] = x;
}

void biot_savart_map(double x0[], double y0[], const double R, const double D, const double theta, const int m, const int n, double *Sreal, double *Simag){
    double B[3],M[3];
    int ind_x,ind_y,counter=0;
	int t;
    
    double c=cos(theta),s=sin(theta);
	
	double RC[Nb_angles],RS[Nb_angles];
	for (t=0;t<Nb_angles;t++){
		RC[t] = R*cos(2*M_PI*t/Nb_angles);
		RS[t] = R*sin(2*M_PI*t/Nb_angles);
	}
    
    for (ind_y=0;ind_y < n; ind_y++){
        for (ind_x=0;ind_x < m; ind_x++){
            M[0] = x0[ind_x];
            M[1] = y0[ind_y];
            M[2] = 0;
            rotation(M,c,s);
            M[1]-=D;
            permutation(M);
            biot_savart(B,M,RC,RS,I);
            inv_permutation(B);
            rotation(B,c,-s);
            Sreal[counter] = B[0];
            Simag[counter] = B[1];
            counter++;
        }
    }
    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *x0,*y0;
    double R,D,theta;
    double *Sreal,*Simag;
    mwSize m,n; 
    
    /* Check for proper number of arguments */
    
    if (nrhs != 5) {
        mexErrMsgTxt("Five inputs required.");
    }
    else if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments."); 
    } 
    
    /* Get the dimensions of x0 and y0. */ 
    
    m = MAX(mxGetM(x0_IN),mxGetN(x0_IN));
    n = MAX(mxGetM(y0_IN),mxGetN(y0_IN));
    
	//mexPrintf("vacuum magnetic permitivity is:  %g\n", mu);
    
    /* Create a matrix for the return argument */ 
    S_OUT = mxCreateDoubleMatrix(m, n, mxCOMPLEX); 
    
    /* Assign pointers to the various parameters */ 
    Sreal = mxGetPr(S_OUT);
    Simag = mxGetPi(S_OUT);
    
    x0    = mxGetPr(x0_IN);
    y0    = mxGetPr(y0_IN);
    R     = *mxGetPr(R_IN);
    D     = *mxGetPr(D_IN);
    theta = *mxGetPr(theta_IN);
        
    /* Do the actual computations in a subroutine */
    biot_savart_map(x0,y0,R,D,theta,m,n,Sreal,Simag);
    
    return;
}
