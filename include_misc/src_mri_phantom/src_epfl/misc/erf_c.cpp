/* This function computes
 *    Erf(R) = 2R/sqrt(pi) * integral from 0 to 1 of exp(-R^2t^2) dt
 * for R real.
 * 
 *  Input: * R: array.
 *
 * This implementation relies no the GSL routines as found in version
 * 1.14.
 *
 * specfunc/erfc.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  J. Theiler (modifications by G. Jungman) */

/*
 * See Hart et al, Computer Approximations, John Wiley and Sons, New York (1968)
 * (This applies only to the erfc8 stuff, which is the part
 *  of the original code that survives. I have replaced much of
 *  the other stuff with Chebyshev fits. These are simpler and
 *  more precise than the original approximations. [GJ])
 *
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

#ifndef M_SQRTPI
#define M_SQRTPI 1.77245385090551588191942755657
#endif

#define	O	plhs[0]
#define	R	prhs[0]

/* Abramowitz+Stegun, 7.1.5 */
static double erfseries(double x)
{
  double coef = x;
  double e    = coef;
  double del;
  int k;
  for (k=1; k<30; ++k) {
    coef *= -x*x/k;
    del   = coef/(2.0*k+1.0);
    e += del;
  }
  return 2.0 / M_SQRTPI * e;
}

struct cheb_series_struct {
  double * c;   /* coefficients                */
  int order;    /* order of expansion          */
  double a;     /* lower interval point        */
  double b;     /* upper interval point        */
  int order_sp; /* effective single precision order */
};
typedef struct cheb_series_struct cheb_series;

/* Chebyshev fit for erfc((t+1)/2), -1 < t < 1
 */
static double erfc_xlt1_data[20] = {
  1.06073416421769980345174155056,
 -0.42582445804381043569204735291,
  0.04955262679620434040357683080,
  0.00449293488768382749558001242,
 -0.00129194104658496953494224761,
 -0.00001836389292149396270416979,
  0.00002211114704099526291538556,
 -5.23337485234257134673693179020e-7,
 -2.78184788833537885382530989578e-7,
  1.41158092748813114560316684249e-8,
  2.72571296330561699984539141865e-9,
 -2.06343904872070629406401492476e-10,
 -2.14273991996785367924201401812e-11,
  2.22990255539358204580285098119e-12,
  1.36250074650698280575807934155e-13,
 -1.95144010922293091898995913038e-14,
 -6.85627169231704599442806370690e-16,
  1.44506492869699938239521607493e-16,
  2.45935306460536488037576200030e-18,
 -9.29599561220523396007359328540e-19
};
static cheb_series erfc_xlt1_cs = {
  erfc_xlt1_data,
  19,
  -1, 1,
  12
};

/* Chebyshev fit for erfc(x) exp(x^2), 1 < x < 5, x = 2t + 3, -1 < t < 1
 */
static double erfc_x15_data[25] = {
  0.44045832024338111077637466616,
 -0.143958836762168335790826895326,
  0.044786499817939267247056666937,
 -0.013343124200271211203618353102,
  0.003824682739750469767692372556,
 -0.001058699227195126547306482530,
  0.000283859419210073742736310108,
 -0.000073906170662206760483959432,
  0.000018725312521489179015872934,
 -4.62530981164919445131297264430e-6,
  1.11558657244432857487884006422e-6,
 -2.63098662650834130067808832725e-7,
  6.07462122724551777372119408710e-8,
 -1.37460865539865444777251011793e-8,
  3.05157051905475145520096717210e-9,
 -6.65174789720310713757307724790e-10,
  1.42483346273207784489792999706e-10,
 -3.00141127395323902092018744545e-11,
  6.22171792645348091472914001250e-12,
 -1.26994639225668496876152836555e-12,
  2.55385883033257575402681845385e-13,
 -5.06258237507038698392265499770e-14,
  9.89705409478327321641264227110e-15,
 -1.90685978789192181051961024995e-15,
  3.50826648032737849245113757340e-16
};
static cheb_series erfc_x15_cs = {
  erfc_x15_data,
  24,
  -1, 1,
  16
};

/* Chebyshev fit for erfc(x) x exp(x^2), 5 < x < 10, x = (5t + 15)/2, -1 < t < 1
 */
static double erfc_x510_data[20] = {
  1.11684990123545698684297865808,
  0.003736240359381998520654927536,
 -0.000916623948045470238763619870,
  0.000199094325044940833965078819,
 -0.000040276384918650072591781859,
  7.76515264697061049477127605790e-6,
 -1.44464794206689070402099225301e-6,
  2.61311930343463958393485241947e-7,
 -4.61833026634844152345304095560e-8,
  8.00253111512943601598732144340e-9,
 -1.36291114862793031395712122089e-9,
  2.28570483090160869607683087722e-10,
 -3.78022521563251805044056974560e-11,
  6.17253683874528285729910462130e-12,
 -9.96019290955316888445830597430e-13,
  1.58953143706980770269506726000e-13,
 -2.51045971047162509999527428316e-14,
  3.92607828989125810013581287560e-15,
 -6.07970619384160374392535453420e-16,
  9.12600607264794717315507477670e-17
};
static cheb_series erfc_x510_cs = {
  erfc_x510_data,
  19,
  -1, 1,
  12
};

static double erfc8_sum(double x)
{
  /* estimates erfc(x) valid for 8 < x < 100 */
  /* This is based on index 5725 in Hart et al */

  static double P[] = {
      2.97886562639399288862,
      7.409740605964741794425,
      6.1602098531096305440906,
      5.019049726784267463450058,
      1.275366644729965952479585264,
      0.5641895835477550741253201704
  };
  static double Q[] = {
      3.3690752069827527677,
      9.608965327192787870698,
      17.08144074746600431571095,
      12.0489519278551290360340491,
      9.396034016235054150430579648,
      2.260528520767326969591866945,
      1.0
  };
  double num=0.0, den=0.0;
  int i;

  num = P[5];
  for (i=4; i>=0; --i) {
      num = x*num + P[i];
  }
  den = Q[6];
  for (i=5; i>=0; --i) {
      den = x*den + Q[i];
  }

  return num/den;
}

inline
static double erfc8(double x)
{
  double e;
  e = erfc8_sum(x);
  e *= exp(-x*x);
  return e;
}

static double
cheb_eval_e(const cheb_series * cs,
            const double x)
{
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  double e = 0.0;

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    e += fabs(y2*temp) + fabs(dd) + fabs(cs->c[j]);
    dd = temp;
  }

  { 
    double temp = d;
    d = y*d - dd + 0.5 * cs->c[0];
    e += fabs(y*temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
  }
  return d;
}

double gsl_sf_erfc_e(double x)
{
  const double ax = fabs(x);
  double e_val, e_err;

  /* CHECK_POINTER(result) */

  if(ax <= 1.0) {
    double t = 2.0*ax - 1.0;
    e_val = cheb_eval_e(&erfc_xlt1_cs, t);
  }
  else if(ax <= 5.0) {
    double ex2 = exp(-x*x);
    double t = 0.5*(ax-3.0);
    e_val = ex2 * cheb_eval_e(&erfc_x15_cs, t);
  }
  else if(ax < 10.0) {
    double exterm = exp(-x*x) / ax;
    double t = (2.0*ax - 15.0)/5.0;
    e_val = exterm * cheb_eval_e(&erfc_x510_cs, t);
  }
  else {
    e_val = erfc8(ax);
  }

  if(x < 0.0) {
    return  2.0 - e_val;
  }
  else {
    return e_val;
  }
}

double gsl_sf_erf_e(double x)
{
  /* CHECK_POINTER(result) */

  if(fabs(x) < 1.0) {
    return erfseries(x);
  }
  else {
    return 1.0 - gsl_sf_erfc_e(x);
  }
}

void myerf_vect(
		double* py,				/* output			*/
		double* px,					/* input			*/
		const unsigned int numel,	/* Number of elements to compute	*/
		const unsigned int inc=1)	/* Increment for the elements (for parallelization) */
{
	/* For each element that must be computed */
    //omp_set_num_threads(4);
#pragma omp parallel for schedule(static,5120)
	for (int i=0;i<numel;i++,px+=inc,py+=inc)
		*py = gsl_sf_erf_e(*px);
	return;
}   /*  myerfzparts_vect */

#ifdef USE_PARALLEL_COMPUTING
struct thread_data{
	double* y;
	double* x;
	unsigned int numel;
};
void *myerf_thread(void *arg)
{
	thread_data *tdata=(thread_data *)arg;
#if PARALLEL_METHOD==1
	myerf_vect(tdata->y, tdata->x, tdata->numel, 1);
#elif PARALLEL_METHOD==2
	myerf_vect(tdata->y, tdata->x, tdata->numel, NUM_THREADS);
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
	mxSetDimensions(O, mxGetDimensions(R), mxGetNumberOfDimensions(R));
	numel=mxGetNumberOfElements(R);
	if (numel==0){return;}
	double* outr = (double*)mxMalloc(numel*sizeof(double));
	double* inR = mxGetPr(R), *inI;
#ifdef USE_PARALLEL_COMPUTING
if (numel < 1200){ /* Threads are inefficient for too few entries */
#endif /* USE_PARALLEL_COMPUTING */
	myerf_vect(outr, inR, numel);
#ifdef USE_PARALLEL_COMPUTING
}else{
	pthread_t       threads[NUM_THREADS];
	thread_data     tdata[NUM_THREADS];
	unsigned int    t;
	unsigned int    numel_per_thread = ceil((float)numel/(float)NUM_THREADS);
	unsigned int    Nspecial = NUM_THREADS-((numel-1)%(numel_per_thread-1))-1;
	double *poutr=outr;
	for(t=0; t<NUM_THREADS; t++){
#if PARALLEL_METHOD==1        
		tdata[t].y   = poutr;
		tdata[t].x      = inR;
#elif PARALLEL_METHOD==2
		tdata[t].y   = (poutr+t);
		tdata[t].x      = (inR+t);
#endif /* PARALLEL_METHOD 2 */
		if (t<NUM_THREADS-Nspecial){
			tdata[t].numel = numel_per_thread;
#if PARALLEL_METHOD==1
			poutr   += numel_per_thread;
			inR     += numel_per_thread;
		}else if(t<NUM_THREADS-1){
			tdata[t].numel = numel_per_thread-1;
			poutr   += numel_per_thread-1;
			inR     += numel_per_thread-1;
#endif /* PARALLEL_METHOD 1 */
		}else{ tdata[t].numel = numel_per_thread-1; }
		pthread_create(&threads[t], NULL, myerf_thread, (void*)&tdata[t]);
	}
	for(t=0; t<NUM_THREADS; t++) {
		pthread_join(threads[t], NULL);
	}
}
#endif /* No USE_PARALLEL_COMPUTING */
mxSetPr(O, outr);
return;
} /*    mexFunction */
