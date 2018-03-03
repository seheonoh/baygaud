#include "baygaud.2dmaps.h"

// 2DBAT user defined functions
// 2D MAPS related



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int save_2dmapsfits(fitsfile *fptr1, char *inputfile, char *outputfile, int nax1, int nax2, velocity_field *twodmap_array)
{
    /* perfect single gaussian */
    int group=1;
    int status=0;
    int nelements=0;
    long fpixel=1;
    FILE *file;
    fitsfile *fptr_out; 

    group=1; status=0;
    nelements = nax1*nax2;

    // remove the output fits if it exists.
    if((file = fopen(outputfile, "r")) != NULL)
    {
        remove(outputfile);
    }

    // open input vf fits
    fits_open_file(&fptr1, inputfile, READONLY, &status);
    // create output fits file
    fits_create_file(&fptr_out, outputfile, &status);
    // copy header
    fits_copy_header(fptr1, fptr_out, &status);
    fits_close_file(fptr_out, &status);

    // write results
    fits_open_file(&fptr_out, outputfile, READWRITE, &status);
    fits_write_img(fptr_out, TFLOAT, fpixel, nelements, twodmap_array[0].data[0], &status);

    fits_close_file(fptr1, &status);
    fits_close_file(fptr_out, &status);
    return 0;
}




int read_2dmaps(TR_ringParameters *TRparam, filename_2dbat *fname, char *card)
{
    int nkeys;
    int status=0;
    int anynul=0;
    int i=0, j=0;

    int vf_e_const_flag=0;
    int mom2_const_flag=0;

    FILE *file_exist;
    fitsfile *fptr1, *fptr2, *fptr3, *fptr4;

    // --------------------------------------------------------------------------------------------- //
    // +++ C. READ VELOCITY FIELDS +++
    // --------------------------------------------------------------------------------------------- //
    // Check if the iput FITS exists...
    if((file_exist = fopen(fname[0].fitsfile_2Dinput_VF, "r")) == NULL)
    {
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        printf("\n%s does not exist in %s\n", fname[0].fitsfile_2Dinput_VF, TRparam[0].wdir);
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        return 0;
    }
    if((file_exist = fopen(fname[0].fitsfile_2Dinput_VF_e, "r")) == NULL)
    {
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        printf("\n%s does not exist in %s\n", fname[0].fitsfile_2Dinput_VF_e, TRparam[0].wdir);
        printf("\nInstead, the user supplied constnat VF_error (km/s) will be used: %.2f\n", TRparam[0].vf_e_user);
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        // copy the input VF name for initializing the VF_e's 2D array below
        vf_e_const_flag = 1; // constant vf_e used
        sprintf(fname[0].fitsfile_2Dinput_VF_e, "%s", fname[0].fitsfile_2Dinput_VF);
    }
    if((file_exist = fopen(fname[0].fitsfile_mom0, "r")) == NULL)
    {
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        printf("\n%s does not exist in %s\n", fname[0].fitsfile_mom0, TRparam[0].wdir);
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        return 0;
    }
    if((file_exist = fopen(fname[0].fitsfile_mom2, "r")) == NULL)
    {
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        printf("\n%s does not exist in %s\n", fname[0].fitsfile_mom2, TRparam[0].wdir);
        printf("\nInstead, the user supplied constnat Vdisp (km/s) will be used: %.2f\n", TRparam[0].vdisp_user);
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        // copy the input VF name for initializing the mom2's 2D array below
        mom2_const_flag = 1; // constant mom2 used
        sprintf(fname[0].fitsfile_mom2, "%s", fname[0].fitsfile_2Dinput_VF);
    }


    int fits_read_key_status=0;
    ffopen(&fptr1, fname[0].fitsfile_2Dinput_VF, READONLY, &status); // open fits
    fits_get_hdrspace(fptr1, &nkeys, NULL, &status);
    // READ NAXIS1
    fits_read_key_status = fits_read_key(fptr1, TINT, "NAXIS1", &TRparam[0].nax1, card, &status);
    fits_read_key_status = fits_read_key(fptr1, TINT, "NAXIS2", &TRparam[0].nax2, card, &status);
    fits_read_key_status = fits_read_key(fptr1, TFLOAT, "CDELT1", &TRparam[0].pixelScale, card, &status);
    if(fits_read_key_status == 202) // if CDELT1 keyword is not in the header, check if CD1_1 is
    {
        status = 0;
        fits_read_key_status = fits_read_key(fptr1, TFLOAT, "CD1_1", &TRparam[0].pixelScale, NULL, &status);
        if(fits_read_key_status == 202) // if CD1_1 keyword is not in the header, put a null value
        {
            TRparam[0].pixelScale = 0.999;
        }
    }
    TRparam[0].pixelScale = fabs(3600*TRparam[0].pixelScale); // in arcsec

    if(status)
    {
        printf("CHECK FITS HEADER: CDELT1 or CD1_1 KEYWORDS DOESN'T EXIST...\n");
        fits_report_error(stderr, status);
        return(status);
    }

    // C-2. Read input velocity field
    ffg2de(fptr1, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF[0].data[0], &anynul, &status); // load VF into array
    ffg2de(fptr1, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF_boxFiltered[0].data[0], &anynul, &status); // load VF into array
    fits_close_file(fptr1, &status); // close fits

    // C-3. Read input sigma field (or mom2)
    ffopen(&fptr2, fname[0].fitsfile_2Dinput_VF_e, READONLY, &status); // open fits
    ffg2de(fptr2, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF_sigma[0].data[0], &anynul, &status); // load VF into array
    ffg2de(fptr2, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF_boxFiltered_sigma[0].data[0], &anynul, &status); // load VF into array
    fits_close_file(fptr2, &status); // close fits

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // if constant vf_e is used (e.g., 1 channel resolution in km/s)
    if(vf_e_const_flag == 1)
    {
        for(i=0; i<TRparam[0].nax1; i++)
        {
            for(j=0; j<TRparam[0].nax2; j++)
            {
                HI_VF_sigma[0].data[j][i] = TRparam[0].vf_e_user; // constant vf_e in km/s
            }   
        }
    }

    // C-4. Read input mom0 field
    ffopen(&fptr3, fname[0].fitsfile_mom0, READONLY, &status); // open fits
    ffg2de(fptr3, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF_mom0[0].data[0], &anynul, &status); // load VF into array
    ffg2de(fptr3, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF_sn[0].data[0], &anynul, &status); // load VF into array
    ffg2de(fptr3, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF_boxFiltered_SN[0].data[0], &anynul, &status); // load VF into array
    fits_close_file(fptr3, &status); // close fits

    // C-5. Read input mom2 field
    ffopen(&fptr4, fname[0].fitsfile_mom2, READONLY, &status); // open fits
    ffg2de(fptr4, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF_mom2[0].data[0], &anynul, &status); // load VF into array
    fits_close_file(fptr4, &status); // close fits

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // if constant vdisp is used (e.g., xxx in km/s)
    if(mom2_const_flag == 1)
    {
        for(i=0; i<TRparam[0].nax1; i++)
        {
            for(j=0; j<TRparam[0].nax2; j++)
            {
                HI_VF_mom2[0].data[j][i] = TRparam[0].vdisp_user; // constant vdisp in km/s
            }   
        }
    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // put a null value to the blank pixels
    for(i=0; i<TRparam[0].nax1; i++)
    {
        for(j=0; j<TRparam[0].nax2; j++)
        {
            // Zero value is often used as a null value in some fits images like SAMI
            if(isnan(HI_VF[0].data[j][i]) || isinf(HI_VF[0].data[j][i]) || fabs(HI_VF[0].data[j][i])<1E-10)
            {
                HI_VF[0].data[j][i] = 1E90;
                HI_VF_boxFiltered[0].data[j][i] = 1E90;
                HI_VF_boxFiltered_SN[0].data[j][i] = 1E90;
                HI_VF_sn[0].data[j][i] = 1E90;
                HI_VF_mom0[0].data[j][i] = 1E90;
                HI_VF_mom2[0].data[j][i] = 1E90;
                HI_VF_sigma[0].data[j][i] = 1E90;
                HI_VF_boxFiltered_sigma[0].data[j][i] = 1E90;
            }
        }
    }

    return 0;
}

int read_3dcube(TR_ringParameters *TRparam, filename_2dbat *fname, char *card)
{
    int nkeys;
    int status=0;
    int anynul=0;
    int i=0, j=0, k=0;

    int vf_e_const_flag=0;
    int mom2_const_flag=0;
	int n_blanks=0;

    FILE *file_exist;
    fitsfile *fptr1, *fptr2, *fptr3, *fptr4;

	const gsl_rng_type * type;
	gsl_rng * r;
	gsl_rng_env_setup();
	type = gsl_rng_default;
	r = gsl_rng_alloc (type);

    // --------------------------------------------------------------------------------------------- //
    // +++ C. READ VELOCITY FIELDS +++
    // --------------------------------------------------------------------------------------------- //
    // Check if the iput FITS exists...
    if((file_exist = fopen(fname[0].raw_cube, "r")) == NULL)
    {
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        printf("\n%s does not exist in %s\n", fname[0].raw_cube, TRparam[0].wdir);
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        return 0;
    }

    int fits_read_key_status=0;
    ffopen(&fptr1, fname[0].raw_cube, READONLY, &status); // open fits
    fits_get_hdrspace(fptr1, &nkeys, NULL, &status);
    // READ NAXIS1
    fits_read_key_status = fits_read_key(fptr1, TINT, "NAXIS1", &TRparam[0].nax1, card, &status);
    fits_read_key_status = fits_read_key(fptr1, TINT, "NAXIS2", &TRparam[0].nax2, card, &status);
    fits_read_key_status = fits_read_key(fptr1, TINT, "NAXIS3", &TRparam[0].nax3, card, &status);
    fits_read_key_status = fits_read_key(fptr1, TFLOAT, "CDELT1", &TRparam[0].pixelScale, card, &status);

    if(fits_read_key_status == 202) // if CDELT1 keyword is not in the header, check if CD1_1 is
    {
        status = 0;
        fits_read_key_status = fits_read_key(fptr1, TFLOAT, "CD1_1", &TRparam[0].pixelScale, NULL, &status);
        if(fits_read_key_status == 202) // if CD1_1 keyword is not in the header, put a null value
        {
            TRparam[0].pixelScale = 0.999;
        }
    }
    TRparam[0].pixelScale = fabs(3600*TRparam[0].pixelScale); // in arcsec

    if(status)
    {
        printf("CHECK FITS HEADER: CDELT1 or CD1_1 KEYWORDS DOESN'T EXIST...\n");
        fits_report_error(stderr, status);
        return(status);
    }

    // C-2. Read input cube

	//printf("%d %d %d\n", TRparam[0].nax1, TRparam[0].nax2, TRparam[0].nax3);
	ffg3de(fptr1, 1, 0, TRparam[0].nax1, TRparam[0].nax2, TRparam[0].nax1, TRparam[0].nax2, TRparam[0].nax3, raw_cube[0].data[0][0], &anynul, &status);

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // put a null value to the blank pixels
/*
    for(i=0; i<TRparam[0].nax1; i++)
    {
        for(j=0; j<TRparam[0].nax2; j++)
        {
			if(isnan(raw_cube[0].data[0][j][i]) || isinf(raw_cube[0].data[0][j][i]) ||
				isnan(raw_cube[0].data[TRparam[0].nax3-1][j][i]) || isinf(raw_cube[0].data[TRparam[0].nax3-1][j][i])) 
			{
				for(k=0; k<TRparam[0].nax3; k++)
				{
					//if(isnan(raw_cube[0].data[k][j][i]) || isinf(raw_cube[0].data[k][j][i])) 
					//{
					//	raw_cube[0].data[k][j][i] = 1E90;
					//}
					//raw_cube_mpi[0].data[i][j][k] = raw_cube[0].data[k][j][i];
				printf("%d %f\n", k, raw_cube[0].data[k][j][i]);
					raw_cube[0].data[k][j][i] = 1E90;
				}
			}
			for(k=0; k<TRparam[0].nax3; k++)
			{
				raw_cube_mpi[0].data[i][j][k] = raw_cube[0].data[k][j][i];
				//printf("%d %f\n", k, raw_cube_mpi[0].data[i][j][k]);
			}
        }
    }
*/

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // put a null value to the blank pixels
	
    for(i=0; i<TRparam[0].nax1; i++)
    {
        for(j=0; j<TRparam[0].nax2; j++)
        {
			n_blanks = 0;
			for(k=0; k<TRparam[0].nax3; k++)
			{
				if(isnan(raw_cube[0].data[k][j][i]) || isinf(raw_cube[0].data[k][j][i]) || fabs(raw_cube[0].data[k][j][i]) < 0.00001)
				{
					raw_cube[0].data[k][j][i] = gsl_ran_gaussian(r, 1E-8);
					n_blanks++;
				}
			}

			if(n_blanks == TRparam[0].nax3) // blank spaxel! put inf value..
			{
				for(k=0; k<TRparam[0].nax3; k++)
				{
					raw_cube[0].data[k][j][i] = 1E90;
				}
			}

			for(k=0; k<TRparam[0].nax3; k++)
			{
				raw_cube_mpi[0].data[i][j][k] = raw_cube[0].data[k][j][i];
				//if(i == 46 && j == 51)
				//printf("%d %d %d %f\n", i, j, k, raw_cube_mpi[0].data[i][j][k]);
			}
        }
    }

	//free(raw_cube[0].data);

    return 0;
}


int read_2d_ref_VF(TR_ringParameters *TRparam, filename_2dbat *fname, char *card)
{
    int nkeys;
    int status=0;
    int anynul=0;
    int i=0, j=0;

    int vf_e_const_flag=0;
    int mom2_const_flag=0;

    FILE *file_exist;
    fitsfile *fptr1, *fptr2, *fptr3, *fptr4;

    // --------------------------------------------------------------------------------------------- //
    // +++ C. READ VELOCITY FIELDS +++
    // --------------------------------------------------------------------------------------------- //
    // Check if the iput FITS exists...
    if((file_exist = fopen(fname[0].fitsfile_trfit_model, "r")) == NULL)
    {
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        printf("\n%s does not exist in %s\n", fname[0].fitsfile_trfit_model, TRparam[0].wdir);
        printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        return 0;
    }


    int fits_read_key_status=0;
    ffopen(&fptr1, fname[0].fitsfile_trfit_model, READONLY, &status); // open fits
    fits_get_hdrspace(fptr1, &nkeys, NULL, &status);
    // READ NAXIS1
    fits_read_key_status = fits_read_key(fptr1, TINT, "NAXIS1", &TRparam[0].nax1, card, &status);
    fits_read_key_status = fits_read_key(fptr1, TINT, "NAXIS2", &TRparam[0].nax2, card, &status);
    //fits_read_key_status = fits_read_key(fptr1, TFLOAT, "CDELT1", &TRparam[0].pixelScale, card, &status);

    if(fits_read_key_status == 202) // if CDELT1 keyword is not in the header, check if CD1_1 is
    {
        status = 0;
        fits_read_key_status = fits_read_key(fptr1, TFLOAT, "CD1_1", &TRparam[0].pixelScale, NULL, &status);
        if(fits_read_key_status == 202) // if CD1_1 keyword is not in the header, put a null value
        {
            TRparam[0].pixelScale = 0.999;
        }
    }
    //TRparam[0].pixelScale = fabs(3600*TRparam[0].pixelScale); // in arcsec

    //if(status)
    //{
    //   printf("CHECK FITS HEADER: CDELT1 or CD1_1 KEYWORDS DOESN'T EXIST...\n");
    //    fits_report_error(stderr, status);
    //    return(status);
    //}

    // C-2. Read input velocity field
    ffg2de(fptr1, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF[0].data[0], &anynul, &status); // load VF into array
    //ffg2de(fptr1, 1, 0, TRparam[0].nax1, TRparam[0].nax1, TRparam[0].nax2, HI_VF_boxFiltered[0].data[0], &anynul, &status); // load VF into array
    fits_close_file(fptr1, &status); // close fits


    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // put a null value to the blank pixels
    for(i=0; i<TRparam[0].nax1; i++)
    {
        for(j=0; j<TRparam[0].nax2; j++)
        {
            // Zero value is often used as a null value in some fits images like SAMI
            if(isnan(HI_VF[0].data[j][i]) || isinf(HI_VF[0].data[j][i]))
            {
                HI_VF[0].data[j][i] = 1E90;
            }
        }
    }
    return 0;
}

// --- End of line

