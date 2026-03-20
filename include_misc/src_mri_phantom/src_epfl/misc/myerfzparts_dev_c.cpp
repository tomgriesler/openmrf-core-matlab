/* This function computes
 *              Erf(R+1j*I) - Erf(R)
 * for R, I reals and R nonzero, with the error function defined as
 *      Erf(z) = 2z/sqrt(pi) * integral from 0 to 1 of exp(-z^2t^2) dt,
 * for any complex z.
 * 
 *  Input: * Z: array
 *              -> if Z is complex, R=real(Z) and I=abs(imag(Z))
 *              -> if Z is real,    R=I=Z .
 *
 * The present implementation is inspired by a matlab implementation from
 * Marcel Leutenegger ( Copyright January 2008, LGPL version 2.1).
 * It uses a series development described in [1]. For details, see [2].
 *
 * Apart the benefits of multi-threading on multi-core units, the main point
 * in this implementation is to evaluate exponentials using multiplications
 * in the summation loops. We also evaluate only the relevant parts.
 * Note, however that:
 * 		- F costs a multiplication when H is evaluated
 * 		- Hr costs a multiplication when Hi is evaluated
 * 		- Gr costs a multiplication when Hi is evaluated
 * 		- Er costs a multiplication when Ei is evaluated
 *
 * [1]  Abramowitz M. and Stegun I.A. (ed.),
 *      "Handbook of Mathematical Functions,"
 *      New York: Dover (1972), Chapter 7, page 299.
 * [2] https://documents.epfl.ch/users/l/le/leuteneg/www/MATLABToolbox/
 *                                                      ErrorFunction.html
 *
 * Contact:     Matthieu Guerquin-Kern,
 *      Biomedical Imaging Group / EPF Lausanne, 31-03-2011 (dd-mm-yyyy)
 */

#define USE_PARALLEL_COMPUTING      /* Implements parallelization using POSIX           */

#include "mex.h"
#include <cmath>
#include <omp.h>

#ifdef USE_PARALLEL_COMPUTING
#include "pthread.h"
#define NUM_THREADS	4
#define PARALLEL_METHOD 1   /* Method 1: (seems more efficient) each thread deals	*
 *				with a distinct block of data			*
 * Method 2: the threads deal with interleaved blocks	*/
#if NUM_THREADS==1
#undef USE_PARALLEL_COMPUTING
#undef PARALLEL_METHOD
#endif /* NUM_THREADS==1 */
#endif /* PARALLEL_COMPUTING */

#define	O	plhs[0]
#define	Z	prhs[0]

#define	MAX(A, B)	((A) > (B) ? (A) : (B))

#ifndef N
#define N 12.17
#define CN 13
#endif /* N */

#if CN<=13
const double nn_4[13]=			/* 0.25*n*n for n from 1 to 13 */
{0.25, 1, 2.25, 4, 6.25, 9, 12.25, 16, 20.25, 25, 30.25, 36, 42.25 };
const double InvExpnn_4[13] =	/*  exp(-0.25*n*n) for n from 1 to 13 */
{   0.77880078307140487848,     0.36787944117144233402,
	0.10539922456186433253,     0.018315638888734178669,
	0.0019304541362277093022,   0.00012340980408667956121,
	4.7851173921290096017e-06,  1.1253517471925911646e-07,
	1.6052280551856116452e-09,  1.3887943864964020896e-11,
	7.287724095819692186e-14,   2.319522830243569634e-16,
	4.4777324417183014957e-19								};
#else /* CN>13 */
#undef CN
#endif /* CN */

#ifndef CN
int CN=ceil(N);
#endif /* No CN */

const double InvSqrtE = sqrt(1/M_E);    /* also exp(-0.5) or 0.60653065971263342426 */

inline int erfz_region(double R2, double aI)
{
	if (aI==0){
		return 0; /* case Z real, nothing to compute: return 0			*/
	}else if (((R2<0.003025)&&(aI>26.63))||((R2>0.003025)&&(aI*aI>712+R2+0.5*log(R2)))){
		return 1; /* G overflows, nothing to compute: return HUGE_VAL	*/
	}else if (R2>36.0){
		return 2; 	/* Compute G										*/
	}else if (R2>33.64){
		if (aI>6.1){
			return 3; /* Compute E and G								*/
		}else{
			return 4; /* Compute E, G, and H							*/
		}
	}else if (R2>0.0){
		if (aI>6.1){
			return 5; /* Compute E, F, and G							*/
		}else{
			return 6; /* Compute E, F, G, and H							*/
		}
	}else { /* R==0 */
		if (aI>6.1){
			return 7; /* Compute Ei and Gi								*/
		}else{
			return 8; /* Compute Ei, Hi, and Gi							*/
		}
	}
}

void compute_G_nmin(
				double n_2,
				double EG,
				double InvA,
				double R2,
				int M,
				double &Gr,
				double &Gi	)
{
	double G;
	for (int counter=1;counter<=M;counter++){	/*	Computation of G: sum for n from nc to nmin    */
		n_2 -= 0.5;					/* n_2	is n/2 and is decreasing (towards nmin/2)					*/
		InvA *= InvSqrtE;			/* InvA is exp(-absI+n/2+1/4)										*
		 	 	 	 	 	 	 	 *		or exp(-absI+(n+1)/2+1/4)*exp(-1/2)							*/
		EG *= InvA;					/* EG	is exp(n*absI-R*R-n*n/4)									*
		 	 	 	 	 	 	 	 *		or exp((n+1)*absI-R*R-(n+1)*(n+1)/4))*exp(-absI+n/2+1/4)	*/
		G  = EG/(n_2*n_2 + R2);		/* G	is exp(n*absI-R2-n*n/4)/(n*n/4 + R*R)						*/
		Gi -= n_2*G;
		Gr += G;
	}
}

void compute_Gi_nmin(
				double n_2,
				double EG,
				double InvA,
				double R2,
				int M,
				double &Gi	)
{
	double G;
	for (int counter=1;counter<=M;counter++){	/*	Computation of G: sum for n from nc to nmin    */
		n_2 -= 0.5;					/* n_2	is n/2 and is decreasing (towards nmin/2)					*/
		InvA *= InvSqrtE;			/* InvA is exp(-absI+n/2+1/4)										*
		 	 	 	 	 	 	 	 *		or exp(-absI+(n+1)/2+1/4)*exp(-1/2)							*/
		EG *= InvA;					/* EG	is exp(n*absI-R*R-n*n/4)									*
		 	 	 	 	 	 	 	 *		or exp((n+1)*absI-R*R-(n+1)*(n+1)/4))*exp(-absI+n/2+1/4)	*/
		G  = EG/(n_2*n_2 + R2);		/* G	is exp(n*absI-R2-n*n/4)/(n*n/4 + R*R)						*/
		Gi -= n_2*G;
	}
}

void compute_F_G_H(
			double nc,
			double EG,
			double a,
			double InvA,
			double R2,
			double InvExpI,
			int M,
			const double *InvExpnn_4,
			const double *nn_4,
			double &F,
			double &Gr,
			double &Gi,
			double &Hr,
			double &Hi	
			)
{
	double n_2 = 0.5*nc, G=EG/(0.25*nc*nc + R2);
	Gi=-n_2*G; /* Here we start with the first term n=nc	*/
	Gr=G;
	compute_G_nmin(n_2,EG,InvA,R2,M,Gr,Gi);
	double H, Hp = 1;
	/* In this loop, G is summed for n from nc+1 to nmax (i.e. nc+CN), H and F are summed for n=1 to CN */
	for (int n=1;n<=CN;n++){    /* for (int n=nc+1;n<=nmax;n++){    */
		H = InvExpnn_4[n-1]/(nn_4[n-1] + R2);	/* H	is exp(-n*n/4)/(n*n/4+R*R)				*/
		F += H;
		Hp *= InvExpI;			/* Hp	is exp(-n*absI)											*/
		H *= Hp;				/* H	is exp(-n*n/4-n*absI)/(n*n/4+R*R);						*/
		Hi += n*H;
		Hr += H;
		n_2 += 0.5;				/* n_2	is n/2 and is increasing (towards nmax/2)				*/
		a *= InvSqrtE;			/* a	is exp(absI-n/2+1/4)									*
		 	 	 	 	 	 	 *		or exp(absI-(n-1)/2+1/4)*exp(-1/2)						*/
		EG *= a;				/* EG	is exp(n*absI-R*R-n*n/4)								*
			 	 	 	 	 	 *		or exp((n-1)*absI-R*R-(n-1)*(n-1)/4)*exp(absI-n/2+1/4)	*/
		G = EG/(n_2*n_2 + R2);	/* G	is exp(n*absI-R*R-n*n/4)/(n*n/4+R*R)					*/
		Gi -= n_2*G;
		Gr += G;
	}
}

void compute_F_G(
			double nc,
			double EG,
			double a,
			double InvA,
			double R2,
			double InvExpI,
			int M,
			const double *InvExpnn_4,
			const double *nn_4,
			double &F,
			double &Gr,
			double &Gi
			)
{
	double n_2 = 0.5*nc, G=EG/(0.25*nc*nc + R2);
	Gi=-n_2*G;
	Gr=G; /* Here we start with the first term n=nc	*/
	compute_G_nmin(n_2,EG,InvA,R2,M,Gr,Gi);
	double H;
	/* In this loop, G is summed for n from nc+1 to nmax (i.e. nc+CN), H and F are summed for n=1 to CN */
	for (int n=1;n<=CN;n++){    /* for (int n=nc+1;n<=nmax;n++){    */
		H = InvExpnn_4[n-1]/(nn_4[n-1] + R2);	/* H	is exp(-n*n/4)/(n*n/4+R*R)				*/
		F += H;
		n_2 += 0.5;				/* n_2	is n/2 and is increasing (towards nmax/2)				*/
		a *= InvSqrtE;			/* a	is exp(absI-n/2+1/4)									*
		 	 	 	 	 	 	 *		or exp(absI-(n-1)/2+1/4)*exp(-1/2)						*/
		EG *= a;				/* EG	is exp(n*absI-R*R-n*n/4)								*
			 	 	 	 	 	 *		or exp((n-1)*absI-R*R-(n-1)*(n-1)/4)*exp(absI-n/2+1/4)	*/
		G = EG/(n_2*n_2 + R2);	/* G	is exp(n*absI-R*R-n*n/4)/(n*n/4+R*R)					*/
		Gi -= n_2*G;
		Gr += G;
	}
}

void compute_G_H(
			double nc,
			double EG,
			double a,
			double InvA,
			double R2,
			double InvExpI,
			int M,
			const double *InvExpnn_4,
			const double *nn_4,
			double &Gr,
			double &Gi,
			double &Hr,
			double &Hi	
			)
{
	double n_2 = 0.5*nc, G=EG/(0.25*nc*nc + R2);
	Gi=-n_2*G;
	Gr=G; /* Here we start with the first term n=nc	*/
	compute_G_nmin(n_2,EG,InvA,R2,M,Gr,Gi);
	double H, Hp = 1;
	/* In this loop, G is summed for n from nc+1 to nmax (i.e. nc+CN), H and F are summed for n=1 to CN */
	for (int n=1;n<=CN;n++){    /* for (int n=nc+1;n<=nmax;n++){    */
		H = InvExpnn_4[n-1]/(nn_4[n-1] + R2);	/* H	is exp(-n*n/4)/(n*n/4+R*R)				*/
		Hp *= InvExpI;			/* Hp	is exp(-n*absI)											*/
		H *= Hp;				/* H	is exp(-n*n/4-n*absI)/(n*n/4+R*R);						*/
		Hi += n*H;
		Hr += H;
		n_2 += 0.5;				/* n_2	is n/2 and is increasing (towards nmax/2)				*/
		a *= InvSqrtE;			/* a	is exp(absI-n/2+1/4)									*
		 	 	 	 	 	 	 *		or exp(absI-(n-1)/2+1/4)*exp(-1/2)						*/
		EG *= a;				/* EG	is exp(n*absI-R*R-n*n/4)								*
			 	 	 	 	 	 *		or exp((n-1)*absI-R*R-(n-1)*(n-1)/4)*exp(absI-n/2+1/4)	*/
		G = EG/(n_2*n_2 + R2);	/* G	is exp(n*absI-R*R-n*n/4)/(n*n/4+R*R)					*/
		Gi -= n_2*G;
		Gr += G;
	}
}

void compute_Gi_Hi(
			double nc,
			double EG,
			double a,
			double InvA,
			double R2,
			double InvExpI,
			int M,
			const double *InvExpnn_4,
			const double *nn_4,
			double &Gi,
			double &Hi	)
{
	double n_2 = 0.5*nc;
	Gi=-n_2*EG/(0.25*nc*nc + R2); /* Here we start with the first term n=nc	*/
	compute_Gi_nmin(n_2,EG,InvA,R2,M,Gi);
	double H, Hp = 1;
	/* In this loop, G is summed for n from nc+1 to nmax (i.e. nc+CN), H and F are summed for n=1 to CN */
	for (int n=1;n<=CN;n++){    /* for (int n=nc+1;n<=nmax;n++){    */
		H = InvExpnn_4[n-1]/(nn_4[n-1] + R2);	/* H	is exp(-n*n/4)/(n*n/4+R*R)				*/
		Hp *= InvExpI;			/* Hp	is exp(-n*absI)											*/
		H *= Hp;				/* H	is exp(-n*n/4-n*absI)/(n*n/4+R*R);						*/
		Hi += n*H;
		n_2 += 0.5;				/* n_2	is n/2 and is increasing (towards nmax/2)				*/
		a *= InvSqrtE;			/* a	is exp(absI-n/2+1/4)									*
		 	 	 	 	 	 	 *		or exp(absI-(n-1)/2+1/4)*exp(-1/2)						*/
		EG *= a;				/* EG	is exp(n*absI-R*R-n*n/4)								*
			 	 	 	 	 	 *		or exp((n-1)*absI-R*R-(n-1)*(n-1)/4)*exp(absI-n/2+1/4)	*/
		Gi -= n_2*EG/(n_2*n_2 + R2);	/* G	is exp(n*absI-R*R-n*n/4)/(n*n/4+R*R)					*/
	}
}

void compute_G(
			double nc,
			double EG,
			double a,
			double InvA,
			double R2,
			double InvExpI,
			int M,
			double &Gr,
			double &Gi
			)
{
	double n_2 = 0.5*nc, G=EG/(0.25*nc*nc + R2);
	Gi=-n_2*G;
	Gr=G; /* Here we start with the first term n=nc	*/
	compute_G_nmin(n_2,EG,InvA,R2,M,Gr,Gi);
	/* In this loop, G is summed for n from nc+1 to nmax (i.e. nc+CN), H and F are summed for n=1 to CN */
	for (int n=1;n<=CN;n++){    /* for (int n=nc+1;n<=nmax;n++){    */
		n_2 += 0.5;				/* n_2	is n/2 and is increasing (towards nmax/2)				*/
		a *= InvSqrtE;			/* a	is exp(absI-n/2+1/4)									*
		 	 	 	 	 	 	 *		or exp(absI-(n-1)/2+1/4)*exp(-1/2)						*/
		EG *= a;				/* EG	is exp(n*absI-R*R-n*n/4)								*
			 	 	 	 	 	 *		or exp((n-1)*absI-R*R-(n-1)*(n-1)/4)*exp(absI-n/2+1/4)	*/
		G = EG/(n_2*n_2 + R2);	/* G	is exp(n*absI-R*R-n*n/4)/(n*n/4+R*R)					*/
		Gi -= n_2*G;
		Gr += G;
	}
}

void compute_Gi(
			double nc,
			double EG,
			double a,
			double InvA,
			double R2,
			double InvExpI,
			int M,
			double &Gi
			)
{
	double n_2 = 0.5*nc;
	Gi=-n_2*EG/(0.25*nc*nc + R2); /* Here we start with the first term n=nc	*/
	compute_Gi_nmin(n_2,EG,InvA,R2,M,Gi);
	/* In this loop, G is summed for n from nc+1 to nmax (i.e. nc+CN), H and F are summed for n=1 to CN */
	for (int n=1;n<=CN;n++){    /* for (int n=nc+1;n<=nmax;n++){    */
		n_2 += 0.5;				/* n_2	is n/2 and is increasing (towards nmax/2)				*/
		a *= InvSqrtE;			/* a	is exp(absI-n/2+1/4)									*
		 	 	 	 	 	 	 *		or exp(absI-(n-1)/2+1/4)*exp(-1/2)						*/
		EG *= a;				/* EG	is exp(n*absI-R*R-n*n/4)								*
			 	 	 	 	 	 *		or exp((n-1)*absI-R*R-(n-1)*(n-1)/4)*exp(absI-n/2+1/4)	*/
		Gi -= n_2*EG/(n_2*n_2 + R2);	/* G	is exp(n*absI-R*R-n*n/4)/(n*n/4+R*R)					*/
	}
}

/* This function computes
 *              Erf(R+1j*I) - Erf(R)
 * for R, I reals and R nonzero, with the error function defined as
 *      Erf(z) = 2z/sqrt(pi) * integral from 0 to 1 of exp(-z^2t^2) dt,
 * for any complex z.
 *
 * The input and output are arrays of doubles.
 * Each element of the array is dealt with separately.
 */
void myerfzparts_scal(
		double &outr,	/* Real part of the output			*/
		double &outi,	/* Imaginary part of the output		*/
		double R,		/* Real part of the input			*/
		double I)		/* Imaginary part of the input		*/
{
#ifndef CN
	/* 	Build lookup tables	*/
	double nn_4[CN], InvExpnn_4[CN];
	for (int n=0;n<CN;n++){ nn_4[n]=(n+1)*(n+1)*0.25;InvExpnn_4[n]=exp(-nn_4[n]);}
#endif /* No CN */
	double R2=R*R, absI=fabs(I);
	int region = 6;//erfz_region(R2,absI);
	if (region==0){ 		/* case Z real			*/
		outr = 0.0;
		outi = 0.0;
	}else if (region==1){ 	/* G overflows but real pa			*/
		outr = HUGE_VAL;
		outi = -HUGE_VAL;
	}else{
		/* First deal with E */		
		double Er, Ei;
		double Cos2RI=cos(2.0*R*absI), Sin2RI=-sin(2.0*R*absI);
		if (region==2){ /* E underflows */
			Ei 	= 0.0;
			Er	= 0.0;		
		}else if (region<7){
			Ei 	= -Sin2RI/R;
			Er	= (1-Cos2RI)/R;	
		}else if(region>6){
			Ei 	= 2*absI;	/* case R=0 */
			Er	= 0.0;	
		}
		double InvExpR2=exp(-R2), InvExpI=exp(-absI);
		int nmin    = MAX(1,floor(2.0*absI-N)), nmax    = ceil(2*absI+N);
		int nc      = MAX(1,nmax-CN),           M       = nc-nmin;
		double A=exp(absI-nc*0.5+0.25), a=A, InvA = 1/(A*InvSqrtE), EGc=exp( nc*absI- R2 -0.25*nc*nc), EG=EGc;
		double Gi=0.0, Gr=0.0, F=0.0, Hr=0.0, Hi=0.0; 
		if (region<4){ /* Compute G alone */
			compute_G(nc, EG, a, InvA, R2, InvExpI, M, Gr, Gi);
		}else if (region==4){	 	/* Compute G and H */
			compute_G_H(nc, EG, a, InvA, R2, InvExpI, M, InvExpnn_4, nn_4, Gr, Gi, Hr, Hi);
		}else if (region==5){	 	/* Compute F and G */
			compute_F_G(nc, EG, a, InvA, R2, InvExpI, M, InvExpnn_4, nn_4, F, Gr, Gi);
		}else if (region==6){	 	/* Compute F, G and H */
			compute_F_G_H(nc, EG, a, InvA, R2, InvExpI, M, InvExpnn_4, nn_4, F, Gr, Gi, Hr, Hi);			
		}else if (region==7){	 	/* Compute Gi */
			compute_Gi(nc, EG, a, InvA, R2, InvExpI, M, Gi);
		}else{ /* region 8 */
			compute_Gi_Hi(nc, EG, a, InvA, R2, InvExpI, M, InvExpnn_4, nn_4, Gi, Hi);
		}
		outr = InvExpR2*(Er + R*F*2.0 - Cos2RI*R*Hr + 0.5*Sin2RI*Hi) - Cos2RI*R*Gr + Sin2RI*Gi;
		outi = InvExpR2*(Ei - Sin2RI*R*Hr - 0.5*Cos2RI*Hi)		     - Cos2RI*Gi + Sin2RI*R*Gr;
	}
	outr *= 0.5*M_1_PI;
	outi *= 0.5*M_1_PI;
	/* we dealt with absI, now we use Erf(conj(Z)) = conj(Erf(Z)) in the case I<0 */
	if (I<0){ outi = -outi;}
	return;
}  /*  myerfzparts_scal */

void myerfzparts_vect(
		double* poutr,				/* Real part of the output			*/
		double* pouti,				/* Imaginary part of the output		*/
		double* pR,					/* Real part of the input			*/
		double* pI,					/* Imaginary part of the input		*/
		const unsigned int numel,	/* Number of elements to compute	*/
		const unsigned int inc=1)	/* Increment for the elements (for parallelization) */
{
	/* For each element that must be computed */
    //omp_set_num_threads(4);
//#pragma omp parallel for schedule(static,5120)
	for (int i=0;i<numel;i++,poutr+=inc,pouti+=inc,pR+=inc,pI+=inc)
		myerfzparts_scal(*poutr,*pouti,*pR,*pI);
		//*poutr = (double)erfz_region(fabs(*pR),fabs(*pI));
		//*pouti = 0.0;
	return;
}   /*  myerfzparts_vect */

#ifdef USE_PARALLEL_COMPUTING
struct thread_data{
	double* outr;
	double* outi;
	double* R;
	double* I;
	unsigned int numel;
};
void *myerfzparts_thread(void *arg)
{
	thread_data *tdata=(thread_data *)arg;
#if PARALLEL_METHOD==1
	myerfzparts_vect(tdata->outr, tdata->outi, tdata->R, tdata->I, tdata->numel, 1);
#elif PARALLEL_METHOD==2
	myerfzparts_vect(tdata->outr, tdata->outi, tdata->R, tdata->I, tdata->numel, NUM_THREADS);
#endif /* PARALLEL_METHOD 2 */
	pthread_exit(NULL);
} /* myerfzparts_thread */
#endif /* USE_PARALLEL_COMPUTING */

/*************** MAIN ****************************/
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	int numel;
	if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");
	if(nrhs=0) mexErrMsgTxt("Wrong number of input arguments.");
	O=mxCreateDoubleMatrix(0, 0, mxREAL);
	mxSetDimensions(O, mxGetDimensions(Z), mxGetNumberOfDimensions(Z));
	numel=mxGetNumberOfElements(Z);
	if (numel==0){return;}
	double* outr = (double*)mxMalloc(numel*sizeof(double));
	double* outi = (double*)mxMalloc(numel*sizeof(double));
	double* inR = mxGetPr(Z), *inI;
	if (mxIsComplex(Z))
	{
		inI = mxGetPi(Z);
	}else{
		inI = inR;        
	}
#ifdef USE_PARALLEL_COMPUTING
if (numel < 1200){ /* Threads are inefficient for too few entries */
#endif /* USE_PARALLEL_COMPUTING */
	myerfzparts_vect(outr, outi, inR, inI, numel);
#ifdef USE_PARALLEL_COMPUTING
}else{
	pthread_t       threads[NUM_THREADS];
	thread_data     tdata[NUM_THREADS];
	unsigned int    t;
	unsigned int    numel_per_thread = ceil((float)numel/(float)NUM_THREADS);
	unsigned int    Nspecial = NUM_THREADS-((numel-1)%(numel_per_thread-1))-1;
	double *poutr=outr, *pouti=outi;
	for(t=0; t<NUM_THREADS; t++){
#if PARALLEL_METHOD==1        
		tdata[t].outr   = poutr;
		tdata[t].outi   = pouti;
		tdata[t].R      = inR;
		tdata[t].I      = inI;
#elif PARALLEL_METHOD==2
		tdata[t].outr   = (poutr+t);
		tdata[t].outi   = (pouti+t);
		tdata[t].R      = (inR+t);
		tdata[t].I      = (inI+t);
#endif /* PARALLEL_METHOD 2 */
		if (t<NUM_THREADS-Nspecial){
			tdata[t].numel = numel_per_thread;
#if PARALLEL_METHOD==1
			poutr   += numel_per_thread;
			pouti   += numel_per_thread;
			inR     += numel_per_thread;
			inI     += numel_per_thread;
		}else if(t<NUM_THREADS-1){
			tdata[t].numel = numel_per_thread-1;
			poutr   += numel_per_thread-1;
			pouti   += numel_per_thread-1;
			inR     += numel_per_thread-1;
			inI     += numel_per_thread-1;
#endif /* PARALLEL_METHOD 1 */
		}else{ tdata[t].numel = numel_per_thread-1; }
		pthread_create(&threads[t], NULL, myerfzparts_thread, (void*)&tdata[t]);
	}
	for(t=0; t<NUM_THREADS; t++) {
		pthread_join(threads[t], NULL);
	}
}
#endif /* No USE_PARALLEL_COMPUTING */
mxSetPi(O, outi);
mxSetPr(O, outr);
return;
} /*    mexFunction */
