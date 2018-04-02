#include "baygaud.memory.h"

// 2DBAT user defined functions
// memory related

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void free_2dmap_arrays(TR_ringParameters *TRparam, multinest_paramters *multinest_param, float *fits_pointer)
{
    free(HI_VF[0].data);
    free(HI_VF_mom0[0].data);
    free(HI_VF_mom2[0].data);
    free(HI_VF_median[0].data);
    free(HI_VF_res[0].data);
    free(HI_VF_sn[0].data);
    free(HI_VF_sn_median[0].data);
    free(HI_VF_sn_res[0].data);
    free(HI_VF_sigma[0].data);
    free(HI_VF_boxFiltered_vlos_ew[0].data);
    free(HI_VF_sigma_res[0].data);
    free(HI_VF_boxFiltered[0].data);
    free(HI_VF_boxFiltered_sigma[0].data);
    free(HI_VF_boxFiltered_sigma_e_norm[0].data);
    free(HI_VF_sigma_geo_flux_weighted[0].data);
    free(HI_VF_boxFiltered_SN[0].data);
    free(HI_VF_boxFiltered_decim0[0].data);
    free(HI_VF_boxFiltered_decim_user[0].data);
    free(HI_VF_geo_radial_angle_w[0].data);
    free(HI_VF_weight_TRfit[0].data);
    free(HI_VF_fract_navail_nall[0].data);
    free(HI_VF_tr_model[0].data);
    free(HI_VF_einasto_halomodel[0].data);
    free(HI_VF_res_input_minus_trfit[0].data);
    free(HI_VF_res_input_minus_einastofit[0].data);
    free(HI_VF_res_trfit_minus_einastofit[0].data);
    free(HI_VF_temp[0].data);
    free(TRparam);
    free(multinest_param);
    free(fits_pointer);
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void malloc_2dmap_arrays(int nax1, int nax2, float *fits_pointer)
{
    int i=0;
    // --------------------------------------------------------------------------------------------- //
    // +++ B. INITIALISATION +++
    // --------------------------------------------------------------------------------------------- //
    // --- B-1. Initialise parameters for velocity field ---    
    // 1. Velocity field area

    // --- A-4. Arrays for 2D velocity fields ---   
    // Dynamic allocation
    HI_VF[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF[0].data[i] = fits_pointer + i * nax2;

    HI_VF_mom0[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_mom0[0].data[i] = fits_pointer + i * nax2;

    HI_VF_mom2[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_mom2[0].data[i] = fits_pointer + i * nax2;

    HI_VF_median[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_median[0].data[i] = fits_pointer + i * nax2;

    HI_VF_res[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_res[0].data[i] = fits_pointer + i * nax2;

    HI_VF_sigma[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_sigma[0].data[i] = fits_pointer + i * nax2;

    HI_VF_boxFiltered_vlos_ew[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_boxFiltered_vlos_ew[0].data[i] = fits_pointer + i * nax2;

    HI_VF_sigma_res[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_sigma_res[0].data[i] = fits_pointer + i * nax2;

    HI_VF_sn[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_sn[0].data[i] = fits_pointer + i * nax2;

    HI_VF_sn_median[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_sn_median[0].data[i] = fits_pointer + i * nax2;

    HI_VF_sn_res[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_sn_res[0].data[i] = fits_pointer + i * nax2;

    HI_VF_geo_radial_angle_w[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_geo_radial_angle_w[0].data[i] = fits_pointer + i * nax2;

    HI_VF_weight_TRfit[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_weight_TRfit[0].data[i] = fits_pointer + i * nax2;

    HI_VF_fract_navail_nall[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_fract_navail_nall[0].data[i] = fits_pointer + i * nax2;

    HI_VF_boxFiltered[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_boxFiltered[0].data[i] = fits_pointer + i * nax2;

    HI_VF_boxFiltered_sigma[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_boxFiltered_sigma[0].data[i] = fits_pointer + i * nax2;

    HI_VF_boxFiltered_sigma_e_norm[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_boxFiltered_sigma_e_norm[0].data[i] = fits_pointer + i * nax2;

    HI_VF_boxFiltered_SN[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_boxFiltered_SN[0].data[i] = fits_pointer + i * nax2;

    HI_VF_boxFiltered_decim0[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_boxFiltered_decim0[0].data[i] = fits_pointer + i * nax2;

    HI_VF_boxFiltered_decim_user[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_boxFiltered_decim_user[0].data[i] = fits_pointer + i * nax2;

    HI_VF_tr_model[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_tr_model[0].data[i] = fits_pointer + i * nax2;

    HI_VF_einasto_halomodel[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_einasto_halomodel[0].data[i] = fits_pointer + i * nax2;

    HI_VF_res_input_minus_trfit[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_res_input_minus_trfit[0].data[i] = fits_pointer + i * nax2;

    HI_VF_res_input_minus_einastofit[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_res_input_minus_einastofit[0].data[i] = fits_pointer + i * nax2;

    HI_VF_res_trfit_minus_einastofit[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_res_trfit_minus_einastofit[0].data[i] = fits_pointer + i * nax2;

    HI_VF_sigma_geo_flux_weighted[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_sigma_geo_flux_weighted[0].data[i] = fits_pointer + i * nax2;

    HI_VF_temp[0].data = malloc(nax1 * sizeof(float *));
    fits_pointer = malloc(nax1 * nax2 * sizeof(float));
    for(i = 0; i < nax1; i++)
        HI_VF_temp[0].data[i] = fits_pointer + i * nax2;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void malloc_2dmap_arrays_bvf(int nax1, int nax2, float *fits_pointer)
{
    int i=0;
	int nfits=0;
    // --------------------------------------------------------------------------------------------- //
    // +++ B. INITIALISATION +++
    // --------------------------------------------------------------------------------------------- //
    // --- B-1. Initialise parameters for velocity field ---    
    // 1. Velocity field area

    // --- A-4. Arrays for 2D velocity fields ---   
    // Dynamic allocation
    HI_VF[0].data = malloc(nax2 * sizeof(float *));
    fits_pointer = malloc(nax2 * nax1 * sizeof(float));
    for(i = 0; i < nax2; i++)
        HI_VF[0].data[i] = fits_pointer + i * nax1;

    // --- A-4. Arrays for 2D velocity fields ---   
    // Dynamic allocation
    HI_VF_boxFiltered[0].data = malloc(nax2 * sizeof(float *));
    fits_pointer = malloc(nax2 * nax1 * sizeof(float));
    for(i = 0; i < nax2; i++)
        HI_VF_boxFiltered[0].data[i] = fits_pointer + i * nax1;

    // --- A-4. Arrays for 2D velocity fields ---   
    // Dynamic allocation
    HI_VF_boxFiltered_sigma[0].data = malloc(nax2 * sizeof(float *));
    fits_pointer = malloc(nax2 * nax1 * sizeof(float));
    for(i = 0; i < nax2; i++)
        HI_VF_boxFiltered_sigma[0].data[i] = fits_pointer + i * nax1;


	// gfit_results!!!!
    // --- A-4. Arrays for 2D velocity fields ---   
    // Dynamic allocation
	for(nfits=0; nfits<120; nfits++)
	{
		BVF[nfits].data = malloc(nax2 * sizeof(float *));
		fits_pointer = malloc(nax2 * nax1 * sizeof(float));
		for(i = 0; i < nax2; i++)
		{
			BVF[nfits].data[i] = fits_pointer + i * nax1;
		}

		BVF_e[nfits].data = malloc(nax2 * sizeof(float *));
		fits_pointer = malloc(nax2 * nax1 * sizeof(float));
		for(i = 0; i < nax2; i++)
		{
			BVF_e[nfits].data[i] = fits_pointer + i * nax1;
		}
	}

	// gfit_results!!!!
    // --- A-4. Arrays for 2D velocity fields ---   
    // Dynamic allocation
	for(nfits=0; nfits<120; nfits++)
	{
		BVF_analysis[nfits].data = malloc(nax2 * sizeof(float *));
		fits_pointer = malloc(nax2 * nax1 * sizeof(float));
		for(i = 0; i < nax2; i++)
		{
			BVF_analysis[nfits].data[i] = fits_pointer + i * nax1;
		}

		BVF_analysis_e[nfits].data = malloc(nax2 * sizeof(float *));
		fits_pointer = malloc(nax2 * nax1 * sizeof(float));
		for(i = 0; i < nax2; i++)
		{
			BVF_analysis_e[nfits].data[i] = fits_pointer + i * nax1;
		}
	}
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void malloc_3dcube_arrays(int nax1, int nax2, int nax3, float *fits_pointer)
{
	int i=0, j=0, k=0;
    float *storage = malloc(nax3 * nax2 * nax1 * sizeof(*storage));
    float *alloc = storage;
    float *storage_mpi = malloc(nax1 * nax2 * nax3 * sizeof(*storage_mpi));
    float *alloc_mpi = storage_mpi;

    raw_cube[0].data = malloc(nax3 * sizeof(*raw_cube[0].data));
    for (k=0; k<nax3; k++)
    {
        raw_cube[0].data[k] = malloc(nax2 * sizeof(**raw_cube[0].data));
        for (j=0; j<nax2; j++)
        {
            raw_cube[0].data[k][j] = alloc;
            alloc += nax1;
        }
    }

    raw_cube_mpi[0].data = malloc(nax1 * sizeof(*raw_cube_mpi[0].data));
    for (i=0; i<nax1; i++)
    {
        raw_cube_mpi[0].data[i] = malloc(nax2 * sizeof(**raw_cube_mpi[0].data));
        for (j=0; j<nax2; j++)
        {
            raw_cube_mpi[0].data[i][j] = alloc_mpi;
            alloc_mpi += nax3;
        }
    }

	//free(storage_mpi);
	//free(storage);
	//free(alloc);
	//free(alloc_mpi);
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void malloc_gfit_params(int nax1, int nax2, int n_gauss, int gfit_params) // 5D array
{
	int i=0, j=0, k=0, l=0;

	// gfit_params
    float *storage_mpi = malloc(nax1 * nax2 * n_gauss * n_gauss * gfit_params * sizeof(*storage_mpi));
    float *alloc_mpi = storage_mpi;
	// gfit_params_e
    float *storage_e_mpi = malloc(nax1 * nax2 * n_gauss * n_gauss * gfit_params * sizeof(*storage_e_mpi));
    float *alloc_e_mpi = storage_e_mpi;

	// allocate gfit_params
    gfit_params_mpi[0].data = malloc(nax1 * sizeof(*gfit_params_mpi[0].data));
    for (i=0; i<nax1; i++)
    {
        gfit_params_mpi[0].data[i] = malloc(nax2 * sizeof(**gfit_params_mpi[0].data));
        for (j=0; j<nax2; j++)
        {
			gfit_params_mpi[0].data[i][j] = malloc(n_gauss * sizeof(***gfit_params_mpi[0].data));
			for (k=0; k<n_gauss; k++)
			{
				gfit_params_mpi[0].data[i][j][k] = malloc(n_gauss * sizeof(****gfit_params_mpi[0].data));
				for (l=0; l<n_gauss; l++)
				{
					gfit_params_mpi[0].data[i][j][k][l] = alloc_mpi;
					alloc_mpi += gfit_params; // gfit_params = 8 (bg(e),A(e),sig(e),X(e) + BIC(e) + STD + S/N)
				}
			}
        }
    }

	// allocate gfit_params_e
    gfit_params_e_mpi[0].data = malloc(nax1 * sizeof(*gfit_params_e_mpi[0].data));
    for (i=0; i<nax1; i++)
    {
        gfit_params_e_mpi[0].data[i] = malloc(nax2 * sizeof(**gfit_params_e_mpi[0].data));
        for (j=0; j<nax2; j++)
        {
			gfit_params_e_mpi[0].data[i][j] = malloc(n_gauss * sizeof(***gfit_params_e_mpi[0].data));
			for (k=0; k<n_gauss; k++)
			{
				gfit_params_e_mpi[0].data[i][j][k] = malloc(n_gauss * sizeof(****gfit_params_e_mpi[0].data));
				for (l=0; l<n_gauss; l++)
				{
					gfit_params_e_mpi[0].data[i][j][k][l] = alloc_e_mpi;
					alloc_e_mpi += gfit_params; // gfit_params = 5 (bg(e),A(e),sig(e),X(e) + BIC(e))
				}
			}
        }
    }

	//free(storage);
	//free(storage_mpi);
	//free(alloc);
	//free(alloc_mpi);
}

void stacksize_()
{
	int res;
	struct rlimit rlim;

	getrlimit(RLIMIT_STACK, &rlim);
	//printf("Before: cur=%d,hard=%d\n",(int)rlim.rlim_cur,(int)rlim.rlim_max);

	rlim.rlim_cur=RLIM_INFINITY;
	rlim.rlim_max=RLIM_INFINITY;
	res=setrlimit(RLIMIT_STACK, &rlim);

	getrlimit(RLIMIT_STACK, &rlim);
	//printf("After: res=%d,cur=%d,hard=%d\n",res,(int)rlim.rlim_cur,(int)rlim.rlim_max);
}


// --- End of line


