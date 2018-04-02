#ifndef __2DBAT_MEMORY_H__
#define __2DBAT_MEMORY_H__

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

#include <mpi.h>

#include "baygaud.cfitsio.h"
#include "baygaud.multinest.h"
#include "baygaud.trfit.h"
#include "baygaud.gsl.h"
#include "baygaud.sort.h"
#include "baygaud.global_params.h"
#include "baygaud.2dmaps.h"
#include "baygaud.etc.h"
#include "baygaud.gfit.h"
#include "baygaud.mpi.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>

// 2DBAT user defined functions
// memory related

void free_2dmap_arrays(TR_ringParameters *TRparam, multinest_paramters *multinest_param, float *fits_pointer);
//
void malloc_2dmap_arrays(int nax1, int nax2, float *fits_pointer);
//
void malloc_2dmap_arrays_bvf(int nax1, int nax2, float *fits_pointer);
//
void malloc_3dcube_arrays(int nax1, int nax2, int nax3, float *fits_pointer);
//
void malloc_gfit_params(int nax1, int nax2, int n_gauss, int gfit_params);
//
void stacksize_();

// --- End of line

#endif

