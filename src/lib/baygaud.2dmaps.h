#ifndef __2DBAT_2DMAPS_H__
#define __2DBAT_2DMAPS_H__

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
#include "baygaud.etc.h"
#include "baygaud.gfit.h"
#include "baygaud.memory.h"
#include "baygaud.mpi.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <float.h>
#include <time.h>

// global 2d arrays : accessible to log_likelihood function in multinest
velocity_field HI_VF[1]; // array for 2D vf : global
velocity_field HI_VF_mom0[1]; // array for 2D mom0 : global
velocity_field HI_VF_mom2[1]; // array for 2D mom2 : global
velocity_field HI_VF_median[1]; // array for 2D vf : global
velocity_field HI_VF_res[1]; // array for 2D vf : global
velocity_field HI_VF_sn[1]; // array for 2D vf : global
velocity_field HI_VF_sn_median[1]; // array for 2D vf : global
velocity_field HI_VF_sn_res[1]; // array for 2D vf : global
velocity_field HI_VF_sigma[1]; // array for 2D vf : global
velocity_field HI_VF_boxFiltered_vlos_ew[1]; // array for 2D vf : global
velocity_field HI_VF_sigma_res[1]; // array for 2D vf : global
velocity_field HI_VF_boxFiltered[1]; // array for 2D vf : global
velocity_field HI_VF_boxFiltered_sigma[1]; // array for 2D vf : global
velocity_field HI_VF_boxFiltered_sigma_e_norm[1]; // array for 2D vf : global
velocity_field HI_VF_boxFiltered_SN[1]; // array for 2D vf : global
velocity_field HI_VF_boxFiltered_decim0[1]; // array for 2D vf : global
velocity_field HI_VF_boxFiltered_decim_user[1]; // array for 2D vf : global
velocity_field HI_VF_geo_radial_angle_w[1]; // array for 2D vf : global
velocity_field HI_VF_weight_TRfit[1]; // array for 2D vf : global
velocity_field HI_VF_fract_navail_nall[1]; // array for 2D vf : global
velocity_field HI_VF_tr_model[1]; // array for 2D vf : global
velocity_field HI_VF_einasto_halomodel[1]; // array for 2D vf : global
velocity_field HI_VF_sigma_geo_flux_weighted[1]; // array for 2D vf : global
velocity_field HI_VF_res_input_minus_trfit[1]; // array for 2D vf : global
velocity_field HI_VF_res_input_minus_einastofit[1]; // array for 2D vf : global
velocity_field HI_VF_res_trfit_minus_einastofit[1]; // array for 2D vf : global
velocity_field HI_VF_temp[1]; // array for 2D vf : global

velocity_field BVF[120]; // array for 2D vf : global
velocity_field BVF_e[120]; // array for 2D vf : global

velocity_field BVF_analysis[120]; // array for 2D vf : global
velocity_field BVF_analysis_e[120]; // array for 2D vf : global

_cube raw_cube[1]; // array for 3D cube : global
_cube raw_cube_mpi[1]; // array for 3D cube : global
_gfit_params gfit_params_mpi[1]; // array for 5D gfit params : global
_gfit_params gfit_params_e_mpi[1]; // array for 5D gfit params : global

// 2DBAT user defined functions
// 2D MAPS related

// save the results in FITS and files
int save_2dmapsfits(fitsfile *fptr1, char *inputfile, char *outputfile, int nax1, int nax2, velocity_field *twodmap_array);
//
int read_2d_ref_VF(TR_ringParameters *TRparam, filename_2dbat *fname, char *card);
//
int read_3dcube(TR_ringParameters *TRparam, filename_2dbat *fname, char *card);

// --- End of line

#endif
