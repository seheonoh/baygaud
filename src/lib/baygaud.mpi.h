#ifndef __2DBAT_MPI_H__
#define __2DBAT_MPI_H__

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
#include "baygaud.memory.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <math.h>
#include <float.h>
#include <time.h>

// 2DBAT user defined functions
// MPI related

// sending or receiving MPI data 
void send_mpi_data(int n_node, int rank, MPI_Request *send_req, MPI_Request *recv_req, TR_ringParameters *TRparam, MPI_Datatype TRparam_mpi,  MPI_Status *mpistatus);
//
void receive_mpiData(int n_node, int rank, MPI_Request *send_req, MPI_Request *recv_req, TR_ringParameters *TRparam, MPI_Datatype TRparam_mpi, velocity_field *HI_VF_boxFiltered, velocity_field *HI_VF_boxFiltered_sigma, velocity_field *HI_VF_boxFiltered_sigma_e_norm, velocity_field *HI_VF_mom2, velocity_field *HI_VF_boxFiltered_decim0, velocity_field *HI_VF, MPI_Status *mpistatus);
//
void receive_mpiData_TRparam(int n_node, int rank, MPI_Request *send_req, MPI_Request *recv_req, TR_ringParameters *TRparam, MPI_Datatype TRparam_mpi, velocity_field *HI_VF_boxFiltered, velocity_field *HI_VF_boxFiltered_sigma, velocity_field *HI_VF_boxFiltered_sigma_e_norm, velocity_field *HI_VF_mom2, velocity_field *HI_VF_boxFiltered_decim0, velocity_field *HI_VF, MPI_Status *mpistatus);
//
MPI_Datatype allocate_mpi_dataset(TR_ringParameters *TRparam, MPI_Datatype *type, int *blocklen);


// --- End of line

#endif
