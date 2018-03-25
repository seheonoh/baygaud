#ifndef __2DBAT_GFIT_H__
#define __2DBAT_GFIT_H__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include <mpi.h>

#include "baygaud.cfitsio.h"
#include "baygaud.multinest.h"
#include "baygaud.trfit.h"
#include "baygaud.gsl.h"
#include "baygaud.sort.h"
#include "baygaud.global_params.h"
#include "baygaud.2dmaps.h"
#include "baygaud.etc.h"
#include "baygaud.memory.h"
#include "baygaud.mpi.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <float.h>
#include <time.h>


// 2DBAT user defined functions
// Gaussian fit related

// Gfit multinest 1D fit
void Gfit_multinest(multinest_paramters *multinest_param, TR_ringParameters *TRparam);
void Gfit_multinest_student(multinest_paramters *multinest_param, TR_ringParameters *TRparam);
void *Gfit_multinest_student_thread(void *thread_args);
void loglikelihood_gfit(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam);
void loglikelihood_gfit_student(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam);
void dumper_Gfits(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, TR_ringParameters *TRparam);
void dumper_Gfits_student(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, TR_ringParameters *TRparam);
double gauss_function(double *Cube, int n_gauss, double x);
double gauss_function_at_a_pixel_gfit_params_mpi(float *****gfit_params, int n_gauss, int _i_, int _j_, double x);
double each_gauss_function_at_a_pixel_gfit_params_mpi(float *****gfit_params, int n_gauss, int nth_gauss, int _i_, int _j_, double x);

void ext_bulk_motions(int xlower, int xupper, int ylower, int yupper, multinest_paramters *multinest_param, TR_ringParameters *TRparam, int sorting_order, float ***HI_cube_mpi, int rank, void *thread_args);
void *ext_bulk_motions_thread(void *thread_args);
int cmp(const void *vp, const void *vq);
int check_thread_status(int pnWorkStatus, int nWaitTime, int rank);
void cleanup(void *arg);

// --- End of line

#endif

