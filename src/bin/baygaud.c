#include "baygaud.main.h"

# define N_PROCESS_MAX 128
# define MAX_TIME_LIMIT 1E6 // SEC

// --- Main program --- //
int main(int argc, char *argv[])
{

	if(argc != 12) // 1:wdir 2:cube.fits (m/s) 3:ref.vf.fits (km/s) 4:xlower 5:ylower 6:xupper 7:yupper
	{
		usage_baygaud();
		exit(0);
	}

    // --------------------------------------------------------------------------------------------- //
    // +++ DECLARE PARAMETRES +++
    // --------------------------------------------------------------------------------------------- //

    // File names for FITS  
    int nkeys, nkeys_2d_refvf;
    char card[FLEN_CARD];
	char dirname[1000];
    int anynul=0;
	struct stat st={0};

	// FITS DATA
	int xlower, ylower, xupper, yupper;
	int sorting_order=1;
	int fits_read_key_status=0;
    FILE *file_exist, *ftemp;
    fitsfile *fptr1, *fptr2, *fptr_2d_refvf;
    float *fits_pointer;

    // Dimensions   
    int i=0, j=0, i0=0, j0=0, k=0, n=0;
	int _ng=0, _ng1=0, _ng2=0, opt_gauss_number=0;


	// THREADS
    int status=0;
	int thread_status;
	int thr_id;
	pthread_t p_thread[N_PROCESS_MAX];
	int output_index;

    clock_t start_time = clock();

    // ETC params
	int mi, mj;

	// MPI
    int size;
    int n_node;
	int tag1=1, tag2=2;
    int rank;

    int n_block;
	int block_length_subcube, stride_subcube;
	int block_length_subcube_gfit_params, stride_subcube_gfit_params;

	int block_length_gaufit, stride_gaufit;
	int block_length_subvf, stride_vf;
	int n_elements, n_stride, x_offset, y_offset;

    int process_X[N_PROCESS_MAX][2]; // process_X[0][0]: xlower, process_X[0][1]: xupper
    int process_Y[N_PROCESS_MAX][2]; // process_Y[0][0]: ylower, process_Y[0][1]: yupper

	double hist_mean_filterbox, hist_std_filterbox;
	float vel_min, vel_max, gauss_model, vel_norm, vel;

	double sn_profile_limit = 0;
	//double channel_resolution = 5200; // in m/s
	double channel_resolution = 0; // in m/s : TBD
	double vesc_limit = 1E9; // in m/s : TBD
	double vel_bulk_limit = 1E9; // in m/s : TBD
	int gauss_number_primary=0, gauss_number_lowest=0, gauss_number_cold=0, gauss_number_warm=0, gauss_number_strong_nonc=0, gauss_number_weak_nonc=0, gauss_number_hvc=0;

	double multinest_time_limit=0;

    MPI_Status mpistatus;
    MPI_Request request[10];
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_node);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    // TR param struct
    TRparam = (TR_ringParameters *)malloc(sizeof (TR_ringParameters) * 1);              // not required by MultiNest, any additional information user wants to pass
    // multinest param struct
    multinest_param = (multinest_paramters *)malloc(sizeof (multinest_paramters) * 1);
    // filename param struct
    fname = (filename_2dbat *)malloc(sizeof (filename_2dbat) * 1);

	//++++++++++++++++++++++++++++++++++++++++++
	// read input from the user
    read_user_input_baygaud(TRparam, multinest_param, fname, argc, argv);
	//++++++++++++++++++++++++++++++++++++++++++
	
	TRparam[0].n_gauss = atoi(argv[4]);	

	/* initialize dynamic 2D array */
	double **vel_sort = (double**)malloc(sizeof(double *) * TRparam[0].n_gauss); // 2D array: TRparam[0].n_gauss x 2
	for(i=0; i<TRparam[0].n_gauss; i++)
	{
	   vel_sort[i] = (double*)malloc(sizeof(double) * 2);
	}

	/* initialize dynamic 2D array */
	double **amp_sort = (double**)malloc(sizeof(double *) * TRparam[0].n_gauss); // 2D array: TRparam[0].n_gauss x 2
	for(i=0; i<TRparam[0].n_gauss; i++)
	{
	   amp_sort[i] = (double*)malloc(sizeof(double) * 2);
	}

	/* initialize dynamic 2D array */
	double **vel_gauss_n_sort = (double**)malloc(sizeof(double *) * TRparam[0].n_gauss); // 2D array: TRparam[0].n_gauss x 2
	for(i=0; i<TRparam[0].n_gauss; i++)
	{
	   vel_gauss_n_sort[i] = (double*)malloc(sizeof(double) * 2);
	}

	/* initialize dynamic 2D array */
	double **disp_sort = (double**)malloc(sizeof(double *) * TRparam[0].n_gauss); // 2D array: TRparam[0].n_gauss x 2
	for(i=0; i<TRparam[0].n_gauss; i++)
	{
	   disp_sort[i] = (double*)malloc(sizeof(double) * 2);
	}

    MPI_Request send_req[4096], recv_req[4096]; // maximum number of cores : 4096
    MPI_Datatype TRparam_mpi;
    TRparam_mpi = allocate_mpi_dataset(TRparam, type, blocklen);    
    MPI_Type_commit(&TRparam_mpi);


	if(rank==0)
	{
		// read 2d input ref velocity field for collecting header info being used for writing fits at the end... 
		ffopen(&fptr_2d_refvf, fname[0].fitsfile_trfit_model, READONLY, &status); // open fits
		fits_get_hdrspace(fptr_2d_refvf, &nkeys_2d_refvf, NULL, &status);

		// read 3d input raw cube fits
		ffopen(&fptr1, fname[0].raw_cube, READONLY, &status); // open fits
		fits_get_hdrspace(fptr1, &nkeys, NULL, &status);
		// READ NAXIS1
		fits_read_key_status = fits_read_key(fptr1, TINT, "NAXIS1", &TRparam[0].nax1, NULL, &status);
		fits_read_key_status = fits_read_key(fptr1, TINT, "NAXIS2", &TRparam[0].nax2, NULL, &status);
		fits_read_key_status = fits_read_key(fptr1, TINT, "NAXIS3", &TRparam[0].nax3, NULL, &status);
		fits_read_key_status = fits_read_key(fptr1, TFLOAT, "CDELT1", &TRparam[0].pixelScale, NULL, &status);
		fits_read_key_status = fits_read_key(fptr1, TFLOAT, "CDELT3", &TRparam[0].Velocity_Seperation, NULL, &status);
		fits_read_key_status = fits_read_key(fptr1, TFLOAT, "CRVAL3", &TRparam[0].Reference_Velocity, NULL, &status);
		fits_read_key_status = fits_read_key(fptr1, TFLOAT, "CRPIX3", &TRparam[0].velocity_offset, NULL, &status);
		TRparam[0].velocity_offset -= 1;
		TRparam[0].N_channels = TRparam[0].nax3;

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
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// SEND TRparam 
    if(rank == 0)
    {
        // --- D-8. Non-blocking sending TRparam_mpi to other processes
        // Wait until all messages have been sent
        // Note that we use requests and statuses starting from index 1
        for (i=1; i<n_node; i++)
        {
            MPI_Isend(TRparam, 1, TRparam_mpi, i, tag1, MPI_COMM_WORLD, &send_req[i]);
            MPI_Wait(&send_req[i], &mpistatus);
        }
	}
    else
    {
        MPI_Irecv(TRparam, 1, TRparam_mpi, 0, tag1, MPI_COMM_WORLD, &recv_req[rank]);
        MPI_Wait(&recv_req[rank], &mpistatus);
	}
	MPI_Barrier(MPI_COMM_WORLD);


/*
            nax1(3)
         ___________
        |__|__|__|__|
nax2(3) |__|__|__|__|  x nax3(10)
        |__|__|__|__|
        |__|__|__|__|

*/
/*
            nax1(2)
         _____
        |__|__|
nlx2(2) |__|__|  x nax3(5)

*/
    // 1. # of blocks (e.g., if we divide the above section into 2 by 2 (i.e., 4 sections) : 4
    // 2. # of elements in each block (4 x nax3 = 40)
    // 3. # of elements between start of each block (if the memory is continuous then, this equals to 4 x nax3 (=40)

	/*++++++++++++++++++++++++++++++++++++++++++*/
	// A. define one sub-section: [xupper-xlower] x [(yupper-ylower)/n_node]
		//xlower = 0;
		//xupper = TRparam[0].nax1;
		//ylower = 0;
		//yupper = TRparam[0].nax2;

	xlower = atoi(argv[5]);
	ylower = atoi(argv[6]);
	xupper = atoi(argv[7]);
	yupper = atoi(argv[8]);
	output_index = atoi(argv[10]);

	double *_filterbox;
	int box_x=xupper-xlower, box_y=yupper-ylower;
	_filterbox = malloc(sizeof(double) * box_x*box_y);

	n_block = (int)(xupper-xlower)/n_node;
	//1. MPI sub cube
	block_length_subcube = (int)(yupper-ylower)*TRparam[0].nax3;
	stride_subcube = TRparam[0].nax2*TRparam[0].nax3;

	//2. MPI sub gaufit
	block_length_gaufit = (int)(yupper-ylower)*4;
	stride_gaufit = TRparam[0].nax2*4;

	//3. MPI sub VF
	block_length_subvf = (int)(yupper-ylower);
	stride_vf = TRparam[0].nax2;

	// MPI data type vector for one sub-section
	MPI_Datatype sub_cube_HIcube;
	MPI_Type_vector(n_block, block_length_subcube, stride_subcube, MPI_FLOAT, &sub_cube_HIcube); // this is for subcube
	MPI_Type_commit(&sub_cube_HIcube);


	//MPI_Isend(&(gfit_params_mpi[0].data[i0][j0][0][0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, i, tag2, MPI_COMM_WORLD, &send_req[i]);


	/*++++++++++++++++++++++++++++++++++++++++++*/
	// B. send one sub-section from the starting address
	// 1. send a portion of the 3D array to rank1: k*j*i
	//  1.buffer-addres, 2.# of sub, 3.sub, 4.rank#, 5.tag, 6.status
	x_offset = (int)(xupper-xlower)/n_node;
	for(i=0; i<n_node; i++)
	{
		// cpu rank
		process_X[i][0] = xlower + x_offset*i;
		process_X[i][1] = xlower + x_offset*(i+1);
		process_Y[i][0] = ylower;
		process_Y[i][1] = yupper;
	}


	// allocate 3d arrays for 3d cube
	malloc_3dcube_arrays(TRparam[0].nax1, TRparam[0].nax2, TRparam[0].nax3, fits_pointer);
	// allocate 5D arrays for saving multiple gfit results
	malloc_gfit_params(TRparam[0].nax1, TRparam[0].nax2, TRparam[0].n_gauss, 10); // 5D array
	// read 2d maps
	malloc_2dmap_arrays_bvf(TRparam[0].nax1, TRparam[0].nax2, fits_pointer);


	// declare MPI data type of gfit_params being used for sending/receiving
	n_block = (int)(xupper-xlower)/n_node;
	//1. MPI sub gfit_params
	block_length_subcube_gfit_params = (int)(yupper-ylower)*TRparam[0].n_gauss*TRparam[0].n_gauss*10;
	stride_subcube_gfit_params = TRparam[0].nax2*TRparam[0].n_gauss*TRparam[0].n_gauss*10;

	// MPI data type vector for one sub-section
	MPI_Datatype sub_gfit_params;
	MPI_Type_vector(n_block, block_length_subcube_gfit_params, stride_subcube_gfit_params, MPI_FLOAT, &sub_gfit_params); // this is for gaufit
	MPI_Type_commit(&sub_gfit_params);


    // ------------------------------------------------------------ //
    // ------------------------------------------------------------ //
    // ------------------------------------------------------------ //
    // ------------------------------------------------------------ //
    if(rank == 0) // master node
    {
		// --------------------------------------------------------------------------------------------- //
		// +++ A. READ INPUT VF KEYWORDS ++W.+
		// --------------------------------------------------------------------------------------------- //

		// read 3d cube
		read_3dcube(TRparam, fname, card);
		// read input ref VF
		read_2d_ref_VF(TRparam, fname, card);

		/*++++++++++++++++++++++++++++++++++++++++++*/
		// send sub-array
		for(n=0; n<n_node-1; n++)
		{
			i0 = process_X[n+1][0];
			j0 = process_Y[n+1][0];
		  // send the input cube and model vf
			MPI_Send(&raw_cube_mpi[0].data[i0][j0][0], 1, sub_cube_HIcube, n+1, 111, MPI_COMM_WORLD);
			MPI_Isend(&(HI_VF[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, n+1, 112, MPI_COMM_WORLD, &send_req[n+1]);
			MPI_Wait(&send_req[n+1], &mpistatus);

			//send_mpi_data(n_node, rank, send_req, recv_req, TRparam, TRparam_mpi, mpistatus);

		  //MPI_Send(&model_VF_mpi[i0][j0], 1, sub_VF, n+1, 111, MPI_COMM_WORLD);
		  //MPI_Send(&perfect_sGfit_mpi[i0][j0][0], 1, sub_cube_gaufit, n+1, 222, MPI_COMM_WORLD);
		}

/*

		int thread_status_master=0; // should be declared
		int _t=0;
		ThreadArgs *thread_args;
		thread_args = (ThreadArgs *)malloc(sizeof(ThreadArgs));
		thread_args->process_X0 = process_X[rank][0];
		thread_args->process_X1 = process_X[rank][1];
		thread_args->process_Y0 = process_Y[rank][0];
		thread_args->process_Y1 = process_Y[rank][1];
		thread_args->status = 0; // not finished
		thread_args->rank = rank; // rank
		thread_args->multinest_time_limit = multinest_time_limit; // sec

		thr_id = pthread_create(&p_thread[rank], NULL, ext_bulk_motions_thread, (void *)thread_args);
		thread_args->status = 1; // not finished

		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// YOU HAVE TO 
		//pthread_detach(p_thread[rank]);
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		while(1)
		{
			_t++;
			sleep(1);
			if(thread_args->status == 1)
			{
				pthread_join(p_thread[rank], (void **)&thread_status_master);
				break;
			}

			if(_t > MAX_TIME_LIMIT)
			{
				printf("foced to exit\n");
				pthread_cleanup_push(cleanup, "thread first handler");
				pthread_cleanup_pop(0);
				pthread_exit((void *)2);
				break;
			}
		}
*/

		ThreadArgs *thread_args;
		thread_args = (ThreadArgs *)malloc(sizeof(ThreadArgs));
		thread_args->multinest_time_limit = atof(argv[11]);; // sec : time limit for each multinest run

		ext_bulk_motions(process_X[rank][0], process_X[rank][1], process_Y[rank][0], process_Y[rank][1], multinest_param, TRparam, sorting_order, raw_cube_mpi[0].data, rank, thread_args);

		/*++++++++++++++++++++++++++++++++++++++++++*/
		// receive gfit_results from sub-nodes
		for(n=0; n<n_node-1; n++)
		{
			i0 = process_X[n+1][0];
			j0 = process_Y[n+1][0];
			MPI_Recv(&gfit_params_mpi[0].data[i0][j0][0][0][0], 1, sub_gfit_params, n+1, 222, MPI_COMM_WORLD, &mpistatus);
		}
		for(n=0; n<n_node-1; n++)
		{
			i0 = process_X[n+1][0];
			j0 = process_Y[n+1][0];
			MPI_Recv(&gfit_params_e_mpi[0].data[i0][j0][0][0][0], 1, sub_gfit_params, n+1, 333, MPI_COMM_WORLD, &mpistatus);
		}

    }
    // ------------------------------------------------------------ //
    else// sub-nodes with ranks
    {
		// read a portion of the 3D array from rank0
		i0 = process_X[rank][0];
		j0 = process_Y[rank][0];
		//printf("rank%d: %d %d\n", rank, i0, j0);
		MPI_Recv(&raw_cube_mpi[0].data[i0][j0][0], 1, sub_cube_HIcube, 0, 111, MPI_COMM_WORLD, &mpistatus);

		MPI_Irecv(&(HI_VF[0].data[0][0]), TRparam[0].nax1*TRparam[0].nax2, MPI_FLOAT, 0, 112, MPI_COMM_WORLD, &recv_req[rank]);
        MPI_Wait(&recv_req[rank], &mpistatus);

/*
		// creat threads for individual ranks
		int thread_status_slave=1; // should be declared
		int _t=0;
		ThreadArgs *thread_args;
		thread_args = (ThreadArgs *)malloc(sizeof(ThreadArgs));
		thread_args->process_X0 = process_X[rank][0];
		thread_args->process_X1 = process_X[rank][1];
		thread_args->process_Y0 = process_Y[rank][0];
		thread_args->process_Y1 = process_Y[rank][1];
		thread_args->status = 0; // not finished
		thread_args->rank = rank; // rank
		thread_args->multinest_time_limit = multinest_time_limit; // sec

		thr_id = pthread_create(&p_thread[rank], NULL, ext_bulk_motions_thread, (void *)thread_args);

		while(1)
		{
			_t++;
			sleep(1);
			if(thread_args->status == 1)
			{
				pthread_join(p_thread[rank], (void **)&thread_status_slave);
				break;
			}

			if(_t > MAX_TIME_LIMIT)
			{
				printf("foced to exit\n");
				pthread_cleanup_push(cleanup, "thread first handler");
				pthread_cleanup_pop(0);
				pthread_exit((void *)2);
				break;
			}
		}
*/

		ThreadArgs *thread_args;
		thread_args = (ThreadArgs *)malloc(sizeof(ThreadArgs));
		thread_args->multinest_time_limit = atof(argv[11]);; // sec : time limit for each multinest run
		ext_bulk_motions(process_X[rank][0], process_X[rank][1], process_Y[rank][0], process_Y[rank][1], multinest_param, TRparam, sorting_order, raw_cube_mpi[0].data, rank, thread_args);

		MPI_Isend(&gfit_params_mpi[0].data[i0][j0][0][0][0], 1, sub_gfit_params, 0, 222, MPI_COMM_WORLD, &send_req[rank]);
		MPI_Isend(&gfit_params_e_mpi[0].data[i0][j0][0][0][0], 1, sub_gfit_params, 0, 333, MPI_COMM_WORLD, &send_req[rank]);
    }

	if(rank == 0)
	{
		for(n=0; n<n_node; n++)
		{
			i0 = process_X[n][0];
			j0 = process_Y[n][0];

			vel_min = (0 - TRparam[0].velocity_offset)*TRparam[0].Velocity_Seperation + TRparam[0].Reference_Velocity;
			vel_max = (TRparam[0].nax3-1 - TRparam[0].velocity_offset)*TRparam[0].Velocity_Seperation + TRparam[0].Reference_Velocity;
		}

		// for re-use of the 2dbat libraries, we use the following arrays for 2D fits files:
		// HI_VF : input reference velocity field
		// BVF_analysis[0] : extracted bulk vf
		// BVF_analysis_e[0] : extracted bulk vf error
		for(i=0; i<TRparam[0].nax1; i++)
		{
			for(j=0; j<TRparam[0].nax2; j++)
			{
				for(k=0; k<90; k++)
				{
					BVF_analysis[k].data[j][i] = 1E90;
					BVF_analysis_e[k].data[j][i] = 1E90;
				}
			}
		}

		//for(i=0; i<TRparam[0].nax1; i++)
		for(i=xlower; i<xupper; i++)
		{
			//for(j=0; j<TRparam[0].nax2; j++)
			for(j=ylower; j<yupper; j++)
			{
				if(!isinf(raw_cube_mpi[0].data[i][j][0]) && !isnan(raw_cube_mpi[0].data[i][j][TRparam[0].nax3-1]) && !isinf(HI_VF[0].data[j][i]) && !isnan(HI_VF[0].data[j][i]))
				{
					opt_gauss_number = gfit_params_mpi[0].data[i][j][0][0][5];

					// +++++++++++++++++++++++++++++++++++++++++++++++++
					// sort vel_sort: observed - reference : better to be re-written for a simplitity: refer to vel_gauss_n_sort below...
					for(_ng=0; _ng<opt_gauss_number; _ng++)
					{
						vel_sort[_ng][0] = 1E9; // default value
						vel_sort[_ng][1] = 1E9; // default large value
					}
					for(_ng=0; _ng<opt_gauss_number; _ng++)
					{
						vel_sort[_ng][0] = gfit_params_mpi[0].data[i][j][opt_gauss_number-1][_ng][3]; // Xo of the current gaussian
						//HI_VF[0].data[j][i] in km/s : input reference velocity field
						vel_sort[_ng][1] = fabs(gfit_params_mpi[0].data[i][j][opt_gauss_number-1][_ng][3] - 1E3*HI_VF[0].data[j][i])*1E3; // absolute velocity difference between Xo and the ref velocity provided: multipled by 1E3 to used the integer comparison function, array2D_comp
					}
					qsort(vel_sort, opt_gauss_number, sizeof(vel_sort[0]), array2D_comp);
					//printf("%d %d ref:%f bulk:%f\n", i, j, HI_VF[0].data[j][i], vel_sort[0][0]);

					// +++++++++++++++++++++++++++++++++++++++++++++++++
					// sort vel_gauss_n_sort: observed - reference : this is for hvc detection.
					for(_ng=0; _ng<opt_gauss_number; _ng++)
					{
						vel_gauss_n_sort[_ng][0] = 1E9; // default value
						vel_gauss_n_sort[_ng][1] = 1E9; // default large value
					}
					for(_ng=0; _ng<opt_gauss_number; _ng++)
					{
						vel_gauss_n_sort[_ng][0] = _ng; // Xo of the current gaussian
						vel_gauss_n_sort[_ng][1] = fabs(gfit_params_mpi[0].data[i][j][opt_gauss_number-1][_ng][3] - 1E3*HI_VF[0].data[j][i])*1E3; // absolute velocity difference between Xo and the ref velocity provided: multipled by 1E3 to used the integer comparison function, array2D_comp
					}
					qsort(vel_gauss_n_sort, opt_gauss_number, sizeof(vel_gauss_n_sort[0]), array2D_comp);
					//printf("%d %d ref:%f bulk:%f\n", i, j, 1E3*HI_VF[0].data[j][i], vel_gauss_n_sort[0][0]);


					// +++++++++++++++++++++++++++++++++++++++++++++++++
					// sort amp: observed
					for(_ng=0; _ng<opt_gauss_number; _ng++)
					{
						amp_sort[_ng][0] = 1E9; // default value
						amp_sort[_ng][1] = 1E9; // default large value
					}
					for(_ng=0; _ng<opt_gauss_number; _ng++)
					{
						amp_sort[_ng][0] = _ng; // gaussian component number
						amp_sort[_ng][1] = gfit_params_mpi[0].data[i][j][opt_gauss_number-1][_ng][1]*1E5; // amplitude: multipled by 1E5 to used the integer comparison function, array2D_comp
					}
					qsort(amp_sort, opt_gauss_number, sizeof(amp_sort[0]), array2D_comp);
					//printf("%d %d ref:%f bulk:%f\n", i, j, HI_VF[0].data[j][i], amp_sort[0][0]);

					// +++++++++++++++++++++++++++++++++++++++++++++++++
					// sort sigma: observed
					for(_ng=0; _ng<opt_gauss_number; _ng++)
					{
						disp_sort[_ng][0] = 1E9; // default value
						disp_sort[_ng][1] = 1E9; // default large value
					}
					for(_ng=0; _ng<opt_gauss_number; _ng++)
					{
						disp_sort[_ng][0] = _ng; // gaussian component number
						disp_sort[_ng][1] = gfit_params_mpi[0].data[i][j][opt_gauss_number-1][_ng][2]*1E3; // dispersion (gaussian sigma term): multipled by 1E3 to used the integer comparison function, array2D_comp
					}
					qsort(disp_sort, opt_gauss_number, sizeof(disp_sort[0]), array2D_comp);
					//printf("%d %d ref:%f bulk:%f\n", i, j, HI_VF[0].data[j][i], disp_sort[0][0]);

					// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
					// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
					// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
					// Extract specific kinematic component
					// ................................................................................................................
					for(_ng=0; _ng<opt_gauss_number; _ng++) // search for the corresponding gaussian component to the found bulk motion
					{

						// *****************************************************************************************
						// 1. bulk velocity field [0 ~ 9] : (1) vel_sort < 0.1 match... it's not limit yet..; (2) sigma > 1 channel; (3) s/n > sn_profile_limit
						// *****************************************************************************************
						if(fabs(gfit_params_mpi[0].data[i][j][opt_gauss_number-1][_ng][3] - vel_sort[0][0]) < 0.1 &&
							gfit_params_mpi[0].data[i][j][opt_gauss_number-1][_ng][2] > channel_resolution &&  // sigma > 1 channel
							gfit_params_mpi[0].data[i][j][opt_gauss_number-1][_ng][7] > sn_profile_limit) // s/n > sn_profile_limit
						{
							for(_ng2=0; _ng2<8; _ng2++)
							{
								// 0:background, 1:amplitude, 2:sigma, 3:centre, 4:bic, 5:none(opt-N-gauss), 6:std, 7:s/n
								BVF_analysis[10*0+_ng2].data[j][i] = gfit_params_mpi[0].data[i][j][opt_gauss_number-1][_ng][_ng2];
								BVF_analysis_e[10*0+_ng2].data[j][i] = gfit_params_e_mpi[0].data[i][j][opt_gauss_number-1][_ng][_ng2];
							}
								BVF_analysis[10*0+2].data[j][i] = BVF_analysis[10*0+2].data[j][i] / 1000.; // in km/s
								BVF_analysis_e[10*0+2].data[j][i] = BVF_analysis_e[10*0+2].data[j][i] / 1000.; // in km/s
								BVF_analysis[10*0+3].data[j][i] = BVF_analysis[10*0+3].data[j][i] / 1000.; // in km/s
								BVF_analysis_e[10*0+3].data[j][i] = BVF_analysis_e[10*0+3].data[j][i] / 1000.; // in km/s
								BVF_analysis[10*0+_ng2+0].data[j][i] = 1E90;
								BVF_analysis[10*0+_ng2+1].data[j][i] = 1E90;
								BVF_analysis_e[10*0+_ng2+0].data[j][i] = 1E90;
								BVF_analysis_e[10*0+_ng2+1].data[j][i] = 1E90;
						}

						if(opt_gauss_number == 1 && // if optimal_n_gauss == 1
							fabs(gfit_params_mpi[0].data[i][j][opt_gauss_number-1][_ng][3]-1E3*HI_VF[0].data[j][i]) < vel_bulk_limit && // in m/s
							gfit_params_mpi[0].data[i][j][opt_gauss_number-1][_ng][2] > channel_resolution &&  // sigma > 1 channel
							gfit_params_mpi[0].data[i][j][opt_gauss_number-1][_ng][7] > sn_profile_limit) // s/n > sn_profile_limit
						{
							for(_ng2=0; _ng2<8; _ng2++)
							{
								// 0:background, 1:amplitude, 2:sigma, 3:centre, 4:bic, 5:none(opt-N-gauss), 6:std, 7:s/n
								BVF_analysis[10*0+_ng2].data[j][i] = gfit_params_mpi[0].data[i][j][opt_gauss_number-1][_ng][_ng2];
								BVF_analysis_e[10*0+_ng2].data[j][i] = gfit_params_e_mpi[0].data[i][j][opt_gauss_number-1][_ng][_ng2];
							}
								BVF_analysis[10*0+2].data[j][i] = BVF_analysis[10*0+2].data[j][i] / 1000.; // in km/s
								BVF_analysis_e[10*0+2].data[j][i] = BVF_analysis_e[10*0+2].data[j][i] / 1000.; // in km/s
								BVF_analysis[10*0+3].data[j][i] = BVF_analysis[10*0+3].data[j][i] / 1000.; // in km/s
								BVF_analysis_e[10*0+3].data[j][i] = BVF_analysis_e[10*0+3].data[j][i] / 1000.; // in km/s
								BVF_analysis[10*0+_ng2+0].data[j][i] = 1E90;
								BVF_analysis[10*0+_ng2+1].data[j][i] = 1E90;
								BVF_analysis_e[10*0+_ng2+0].data[j][i] = 1E90;
								BVF_analysis_e[10*0+_ng2+1].data[j][i] = 1E90;
						}

						// *****************************************************************************************
						// 2. perfect single gaussian fit [10 ~ 19] : (1) optimal_n_gauss == 1; (2) sigma > 1 channel; (3) s/n > sn_profile_limit
						// *****************************************************************************************
						if(opt_gauss_number == 1 &&
							gfit_params_mpi[0].data[i][j][opt_gauss_number-1][_ng][2] > channel_resolution &&  // sigma > 1 channel
							gfit_params_mpi[0].data[i][j][opt_gauss_number-1][_ng][7] > sn_profile_limit) // s/n > sn_profile_limit
						{
							for(_ng2=0; _ng2<8; _ng2++)
							{
								// 0:background, 1:amplitude, 2:sigma, 3:centre, 4:bic, 5:none(opt-N-gauss), 6:std, 7:s/n
								BVF_analysis[10*1+_ng2].data[j][i] = gfit_params_mpi[0].data[i][j][opt_gauss_number-1][_ng][_ng2];
								BVF_analysis_e[10*1+_ng2].data[j][i] = gfit_params_e_mpi[0].data[i][j][opt_gauss_number-1][_ng][_ng2];
							}
								BVF_analysis[10*1+2].data[j][i] = BVF_analysis[10*1+2].data[j][i] / 1000.; // in km/s
								BVF_analysis_e[10*1+2].data[j][i] = BVF_analysis_e[10*1+2].data[j][i] / 1000.; // in km/s
								BVF_analysis[10*1+3].data[j][i] = BVF_analysis[10*1+3].data[j][i] / 1000.; // in km/s
								BVF_analysis_e[10*1+3].data[j][i] = BVF_analysis_e[10*1+3].data[j][i] / 1000.; // in km/s
								BVF_analysis[10*1+_ng2+0].data[j][i] = 1E90;
								BVF_analysis[10*1+_ng2+1].data[j][i] = 1E90;
								BVF_analysis_e[10*1+_ng2+0].data[j][i] = 1E90;
								BVF_analysis_e[10*1+_ng2+1].data[j][i] = 1E90;
						}

						// *****************************************************************************************
						// 3. default single gaussian fit [20 ~ 29] : (1) optimal_n_gauss > 1; (2) sigma > 1 channel; (3) s/n > sn_profile_limit
						// *****************************************************************************************
						if(gfit_params_mpi[0].data[i][j][0][0][2] > channel_resolution &&  // sigma > 1 channel
							gfit_params_mpi[0].data[i][j][0][0][7] > sn_profile_limit) // s/n > sn_profile_limit
						{
							for(_ng2=0; _ng2<8; _ng2++)
							{
								// 0:background, 1:amplitude, 2:sigma, 3:centre, 4:bic, 5:none(opt-N-gauss), 6:std, 7:s/n
								BVF_analysis[10*2+_ng2].data[j][i] = gfit_params_mpi[0].data[i][j][0][0][_ng2];
								BVF_analysis_e[10*2+_ng2].data[j][i] = gfit_params_e_mpi[0].data[i][j][0][0][_ng2];
							}
								BVF_analysis[10*2+2].data[j][i] = BVF_analysis[10*2+2].data[j][i] / 1000.; // in km/s
								BVF_analysis_e[10*2+2].data[j][i] = BVF_analysis_e[10*2+2].data[j][i] / 1000.; // in km/s
								BVF_analysis[10*2+3].data[j][i] = BVF_analysis[10*2+3].data[j][i] / 1000.; // in km/s
								BVF_analysis_e[10*2+3].data[j][i] = BVF_analysis_e[10*2+3].data[j][i] / 1000.; // in km/s
								BVF_analysis[10*2+_ng2+0].data[j][i] = 1E90;
								BVF_analysis[10*2+_ng2+1].data[j][i] = 1E90;
								BVF_analysis_e[10*2+_ng2+0].data[j][i] = 1E90;
								BVF_analysis_e[10*2+_ng2+1].data[j][i] = 1E90;
						}

						// *****************************************************************************************
						// 4. primary gaussian component [30 ~ 39] : (1) ; (2) sigma > 1 channel; (3) s/n > sn_profile_limit ; (4) amp is the largest
						// *****************************************************************************************
						gauss_number_primary = amp_sort[opt_gauss_number-1][0]; // sorted as above // gaussian component number whose amp is the biggest among the optimally decomposed ones...
						if(opt_gauss_number > 1 && // at least the number of gaussian is larger than 1
							gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_primary][2] > channel_resolution &&  // sigma > 1 channel
							gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_primary][7] > sn_profile_limit) // s/n > sn_profile_limit
						{
							for(_ng2=0; _ng2<8; _ng2++)
							{
								// 0:background, 1:amplitude, 2:sigma, 3:centre, 4:bic, 5:none(opt-N-gauss), 6:std, 7:s/n
								BVF_analysis[10*3+_ng2].data[j][i] = gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_primary][_ng2];
								BVF_analysis_e[10*3+_ng2].data[j][i] = gfit_params_e_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_primary][_ng2];
							}
								BVF_analysis[10*3+2].data[j][i] = BVF_analysis[10*3+2].data[j][i] / 1000.; // in km/s
								BVF_analysis_e[10*3+2].data[j][i] = BVF_analysis_e[10*3+2].data[j][i] / 1000.; // in km/s
								BVF_analysis[10*3+3].data[j][i] = BVF_analysis[10*3+3].data[j][i] / 1000.; // in km/s
								BVF_analysis_e[10*3+3].data[j][i] = BVF_analysis_e[10*3+3].data[j][i] / 1000.; // in km/s
								BVF_analysis[10*3+_ng2+0].data[j][i] = 1E90;
								BVF_analysis[10*3+_ng2+1].data[j][i] = 1E90;
								BVF_analysis_e[10*3+_ng2+0].data[j][i] = 1E90;
								BVF_analysis_e[10*3+_ng2+1].data[j][i] = 1E90;
						}

						// *****************************************************************************************
						// 5. cold gaussian component [40 ~ 49] : (1) opt_gauss_number > 1; (2) sigma1 & sigma2 > 1 channel + sigma1 < sigma2; (3) s/n > sn_profile_limit;
						// *****************************************************************************************
						gauss_number_cold = disp_sort[0][0]; // sorted as above // gaussian component number whose dispersion is the smallest among the optimally decomposed ones...
						if(opt_gauss_number > 1 && // at least the number of gaussian is larger than 1
							gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_cold][2] > channel_resolution &&  // sigma > 1 channel
							gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_cold][7] > sn_profile_limit) // s/n > sn_profile_limit
						{
							for(_ng2=0; _ng2<8; _ng2++)
							{
								// 0:background, 1:amplitude, 2:sigma, 3:centre, 4:bic, 5:none(opt-N-gauss), 6:std, 7:s/n
								BVF_analysis[10*4+_ng2].data[j][i] = gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_cold][_ng2];
								BVF_analysis_e[10*4+_ng2].data[j][i] = gfit_params_e_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_cold][_ng2];
							}
								BVF_analysis[10*4+2].data[j][i] = BVF_analysis[10*4+2].data[j][i] / 1000.; // in km/s
								BVF_analysis_e[10*4+2].data[j][i] = BVF_analysis_e[10*4+2].data[j][i] / 1000.; // in km/s
								BVF_analysis[10*4+3].data[j][i] = BVF_analysis[10*4+3].data[j][i] / 1000.; // in km/s
								BVF_analysis_e[10*4+3].data[j][i] = BVF_analysis_e[10*4+3].data[j][i] / 1000.; // in km/s
								BVF_analysis[10*4+_ng2+0].data[j][i] = 1E90;
								BVF_analysis[10*4+_ng2+1].data[j][i] = 1E90;
								BVF_analysis_e[10*4+_ng2+0].data[j][i] = 1E90;
								BVF_analysis_e[10*4+_ng2+1].data[j][i] = 1E90;
						}

						// *****************************************************************************************
						// 6. warm gaussian component [50 ~ 59] : (1) opt_gauss_number > 1; (2) sigma1 & sigma2 > 1 channel + sigma1 > sigma2; (3) s/n > sn_profile_limit;
						// *****************************************************************************************
						gauss_number_warm = disp_sort[opt_gauss_number-1][0]; // sorted as above // gaussian component number whose dispersion is the largest among the optimally decomposed ones...
						if(opt_gauss_number > 1 && // at least the number of gaussian is larger than 1
							gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_warm][2] > channel_resolution &&  // sigma > 1 channel
							gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_warm][7] > sn_profile_limit) // s/n > sn_profile_limit
						{
							for(_ng2=0; _ng2<8; _ng2++)
							{
								// 0:background, 1:amplitude, 2:sigma, 3:centre, 4:bic, 5:none(opt-N-gauss), 6:std, 7:s/n
								BVF_analysis[10*5+_ng2].data[j][i] = gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_warm][_ng2];
								BVF_analysis_e[10*5+_ng2].data[j][i] = gfit_params_e_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_warm][_ng2];
							}
								BVF_analysis[10*5+2].data[j][i] = BVF_analysis[10*5+2].data[j][i] / 1000.; // in km/s
								BVF_analysis_e[10*5+2].data[j][i] = BVF_analysis_e[10*5+2].data[j][i] / 1000.; // in km/s
								BVF_analysis[10*5+3].data[j][i] = BVF_analysis[10*5+3].data[j][i] / 1000.; // in km/s
								BVF_analysis_e[10*5+3].data[j][i] = BVF_analysis_e[10*5+3].data[j][i] / 1000.; // in km/s
								BVF_analysis[10*5+_ng2+0].data[j][i] = 1E90;
								BVF_analysis[10*5+_ng2+1].data[j][i] = 1E90;
								BVF_analysis_e[10*5+_ng2+0].data[j][i] = 1E90;
								BVF_analysis_e[10*5+_ng2+1].data[j][i] = 1E90;
						}

						// *****************************************************************************************
						// 7. strong nonc. gaussian component [60 ~ 69] : (1) opt_gauss_number > 1; (2) not the bulk but the primary; (3) sigma > 1 channel; (4) s/n > sn_profile_limit;
						// *****************************************************************************************
						// first, find the primary
						gauss_number_primary = amp_sort[opt_gauss_number-1][0]; // sorted as above // gaussian component number whose amp is the biggest among the optimally decomposed ones...
						// second, check it the primary is bulk motion whether or not...
						if(fabs(gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_primary][3] - vel_sort[0][0]) > 0.1) // if not the bulk motion
						{
							gauss_number_strong_nonc = gauss_number_primary; // sorted as above // gaussian component number whose velocitiy deviation is the largest from the bulk among the optimally decomposed ones...
							if(opt_gauss_number > 1 && // at least the number of gaussian is larger than 1
								gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_strong_nonc][2] > channel_resolution &&  // sigma > 1 channel
								gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_strong_nonc][7] > sn_profile_limit) // s/n > sn_profile_limit
							{
								for(_ng2=0; _ng2<8; _ng2++)
								{
									// 0:background, 1:amplitude, 2:sigma, 3:centre, 4:bic, 5:none(opt-N-gauss), 6:std, 7:s/n
									BVF_analysis[10*6+_ng2].data[j][i] = gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_strong_nonc][_ng2];
									BVF_analysis_e[10*6+_ng2].data[j][i] = gfit_params_e_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_strong_nonc][_ng2];
								}
									BVF_analysis[10*6+2].data[j][i] = BVF_analysis[10*6+2].data[j][i] / 1000.; // in km/s
									BVF_analysis_e[10*6+2].data[j][i] = BVF_analysis_e[10*6+2].data[j][i] / 1000.; // in km/s
									BVF_analysis[10*6+3].data[j][i] = BVF_analysis[10*6+3].data[j][i] / 1000.; // in km/s
									BVF_analysis_e[10*6+3].data[j][i] = BVF_analysis_e[10*6+3].data[j][i] / 1000.; // in km/s
									BVF_analysis[10*6+_ng2+0].data[j][i] = 1E90;
									BVF_analysis[10*6+_ng2+1].data[j][i] = 1E90;
									BVF_analysis_e[10*6+_ng2+0].data[j][i] = 1E90;
									BVF_analysis_e[10*6+_ng2+1].data[j][i] = 1E90;
							}

							/*
							if(opt_gauss_number == 1 && // if optimal_n_gauss == 1
								fabs(gfit_params_mpi[0].data[i][j][opt_gauss_number-1][0][3]-1E3*HI_VF[0].data[j][i]) > vel_bulk_limit && // in m/s
								gfit_params_mpi[0].data[i][j][opt_gauss_number-1][0][2] > channel_resolution &&  // sigma > 1 channel
								gfit_params_mpi[0].data[i][j][opt_gauss_number-1][0][7] > sn_profile_limit) // s/n > sn_profile_limit
							{
								for(_ng2=0; _ng2<8; _ng2++)
								{
									// 0:background, 1:amplitude, 2:sigma, 3:centre, 4:bic, 5:none(opt-N-gauss), 6:std, 7:s/n
									BVF_analysis[10*6+_ng2].data[j][i] = gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_strong_nonc][_ng2];
									BVF_analysis_e[10*6+_ng2].data[j][i] = gfit_params_e_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_strong_nonc][_ng2];
								}
									BVF_analysis[10*6+2].data[j][i] = BVF_analysis[10*6+2].data[j][i] / 1000.; // in km/s
									BVF_analysis_e[10*6+2].data[j][i] = BVF_analysis_e[10*6+2].data[j][i] / 1000.; // in km/s
									BVF_analysis[10*6+3].data[j][i] = BVF_analysis[10*6+3].data[j][i] / 1000.; // in km/s
									BVF_analysis_e[10*6+3].data[j][i] = BVF_analysis_e[10*6+3].data[j][i] / 1000.; // in km/s
									BVF_analysis[10*6+_ng2+0].data[j][i] = 1E90;
									BVF_analysis[10*6+_ng2+1].data[j][i] = 1E90;
									BVF_analysis_e[10*6+_ng2+0].data[j][i] = 1E90;
									BVF_analysis_e[10*6+_ng2+1].data[j][i] = 1E90;
							}
							*/
						}

						// *****************************************************************************************
						// 8. weak nonc. gaussian component [70 ~ 79] : (1) opt_gauss_number > 1; (2) not the bulk but the secondary (3) sigma > 1 channel; (4) s/n > sn_profile_limit;
						// *****************************************************************************************
						// first, find the primary
						gauss_number_lowest = amp_sort[0][0]; // sorted as above // gaussian component number whose amp is the lowest among the optimally decomposed ones...
						// second, check it the lowest is bulk motion whether or not...
						if(fabs(gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_lowest][3] - vel_sort[0][0]) > 0.1) // if not the bulk motion
						{
							gauss_number_weak_nonc = gauss_number_lowest; // sorted as above // gaussian component number whose velocitiy deviation is the largest from the bulk among the optimally decomposed ones...
							if(opt_gauss_number > 1 && // at least the number of gaussian is larger than 1
								gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_weak_nonc][2] > channel_resolution &&  // sigma > 1 channel
								gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_weak_nonc][7] > sn_profile_limit) // s/n > sn_profile_limit
							{
								for(_ng2=0; _ng2<8; _ng2++)
								{
									// 0:background, 1:amplitude, 2:sigma, 3:centre, 4:bic, 5:none(opt-N-gauss), 6:std, 7:s/n
									BVF_analysis[10*7+_ng2].data[j][i] = gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_weak_nonc][_ng2];
									BVF_analysis_e[10*7+_ng2].data[j][i] = gfit_params_e_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_weak_nonc][_ng2];
								}
									BVF_analysis[10*7+2].data[j][i] = BVF_analysis[10*7+2].data[j][i] / 1000.; // in km/s
									BVF_analysis_e[10*7+2].data[j][i] = BVF_analysis_e[10*7+2].data[j][i] / 1000.; // in km/s
									BVF_analysis[10*7+3].data[j][i] = BVF_analysis[10*7+3].data[j][i] / 1000.; // in km/s
									BVF_analysis_e[10*7+3].data[j][i] = BVF_analysis_e[10*7+3].data[j][i] / 1000.; // in km/s
									BVF_analysis[10*7+_ng2+0].data[j][i] = 1E90;
									BVF_analysis[10*7+_ng2+1].data[j][i] = 1E90;
									BVF_analysis_e[10*7+_ng2+0].data[j][i] = 1E90;
									BVF_analysis_e[10*7+_ng2+1].data[j][i] = 1E90;
							}

							/*
							if(opt_gauss_number == 1 && // if optimal_n_gauss == 1
								fabs(gfit_params_mpi[0].data[i][j][opt_gauss_number-1][0][3]-1E3*HI_VF[0].data[j][i]) > vel_bulk_limit && // in m/s
								gfit_params_mpi[0].data[i][j][opt_gauss_number-1][0][2] > channel_resolution &&  // sigma > 1 channel
								gfit_params_mpi[0].data[i][j][opt_gauss_number-1][0][7] > sn_profile_limit) // s/n > sn_profile_limit
							{
								for(_ng2=0; _ng2<8; _ng2++)
								{
									// 0:background, 1:amplitude, 2:sigma, 3:centre, 4:bic, 5:none(opt-N-gauss), 6:std, 7:s/n
									BVF_analysis[10*7+_ng2].data[j][i] = gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_weak_nonc][_ng2];
									BVF_analysis_e[10*7+_ng2].data[j][i] = gfit_params_e_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_weak_nonc][_ng2];
								}
									BVF_analysis[10*7+2].data[j][i] = BVF_analysis[10*7+2].data[j][i] / 1000.; // in km/s
									BVF_analysis_e[10*7+2].data[j][i] = BVF_analysis_e[10*7+2].data[j][i] / 1000.; // in km/s
									BVF_analysis[10*7+3].data[j][i] = BVF_analysis[10*7+3].data[j][i] / 1000.; // in km/s
									BVF_analysis_e[10*7+3].data[j][i] = BVF_analysis_e[10*7+3].data[j][i] / 1000.; // in km/s
									BVF_analysis[10*7+_ng2+0].data[j][i] = 1E90;
									BVF_analysis[10*7+_ng2+1].data[j][i] = 1E90;
									BVF_analysis_e[10*7+_ng2+0].data[j][i] = 1E90;
									BVF_analysis_e[10*7+_ng2+1].data[j][i] = 1E90;
							}
							*/
						}

						// *****************************************************************************************
						// 9. HVC gaussian component [80 ~ 89] : (1) del_v > vesc_limit; (2) sigma > 1 channel; (3) s/n > sn_profile_limit;
						// *****************************************************************************************
						if(opt_gauss_number == 1)
						{
							gauss_number_hvc = 0;
							if(fabs(gfit_params_mpi[0].data[i][j][0][0][3] - 1E3*HI_VF[0].data[j][i]) > vesc_limit) // if not the bulk motion
							{
								if(gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_hvc][2] > channel_resolution &&  // sigma > 1 channel
									gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_hvc][7] > sn_profile_limit) // s/n > sn_profile_limit
								{
									for(_ng2=0; _ng2<8; _ng2++)
									{
										// 0:background, 1:amplitude, 2:sigma, 3:centre, 4:bic, 5:none(opt-N-gauss), 6:std, 7:s/n
										BVF_analysis[10*8+_ng2].data[j][i] = gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_hvc][_ng2];
										BVF_analysis_e[10*8+_ng2].data[j][i] = gfit_params_e_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_hvc][_ng2];
									}
										BVF_analysis[10*8+2].data[j][i] = BVF_analysis[10*8+2].data[j][i] / 1000.; // in km/s
										BVF_analysis_e[10*8+2].data[j][i] = BVF_analysis_e[10*8+2].data[j][i] / 1000.; // in km/s
										BVF_analysis[10*8+3].data[j][i] = BVF_analysis[10*8+3].data[j][i] / 1000.; // in km/s
										BVF_analysis_e[10*8+3].data[j][i] = BVF_analysis_e[10*8+3].data[j][i] / 1000.; // in km/s
										BVF_analysis[10*8+_ng2+0].data[j][i] = 1E90;
										BVF_analysis[10*8+_ng2+1].data[j][i] = 1E90;
										BVF_analysis_e[10*8+_ng2+0].data[j][i] = 1E90;
										BVF_analysis_e[10*8+_ng2+1].data[j][i] = 1E90;
								}
							}
						}
						else
						{
							gauss_number_hvc = vel_gauss_n_sort[opt_gauss_number-1][0]; // sorted as above // gaussian component number whose amp is the biggest among the optimally decomposed ones...
							if(fabs(gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_hvc][3] - 1E3*HI_VF[0].data[j][i]) > vesc_limit) // if not the bulk motion
							{
								if(gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_hvc][2] > channel_resolution &&  // sigma > 1 channel
									gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_hvc][7] > sn_profile_limit) // s/n > sn_profile_limit
								{
									for(_ng2=0; _ng2<8; _ng2++)
									{
										// 0:background, 1:amplitude, 2:sigma, 3:centre, 4:bic, 5:none(opt-N-gauss), 6:std, 7:s/n
										BVF_analysis[10*8+_ng2].data[j][i] = gfit_params_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_hvc][_ng2];
										BVF_analysis_e[10*8+_ng2].data[j][i] = gfit_params_e_mpi[0].data[i][j][opt_gauss_number-1][gauss_number_hvc][_ng2];
									}
										BVF_analysis[10*8+2].data[j][i] = BVF_analysis[10*8+2].data[j][i] / 1000.; // in km/s
										BVF_analysis_e[10*8+2].data[j][i] = BVF_analysis_e[10*8+2].data[j][i] / 1000.; // in km/s
										BVF_analysis[10*8+3].data[j][i] = BVF_analysis[10*8+3].data[j][i] / 1000.; // in km/s
										BVF_analysis_e[10*8+3].data[j][i] = BVF_analysis_e[10*8+3].data[j][i] / 1000.; // in km/s
										BVF_analysis[10*8+_ng2+0].data[j][i] = 1E90;
										BVF_analysis[10*8+_ng2+1].data[j][i] = 1E90;
										BVF_analysis_e[10*8+_ng2+0].data[j][i] = 1E90;
										BVF_analysis_e[10*8+_ng2+1].data[j][i] = 1E90;
								}
							}
						}

					}

					// ++++++++++++++++++++++++++++++++++++++++++++++++++++++
					// ++++++++++++++++++++++++++++++++++++++++++++++++++++++
					// ++++++++++++++++++++++++++++++++++++++++++++++++++++++
					// slicing all the fit results from multiple gaussian fits
					for(_ng=0; _ng<opt_gauss_number; _ng++)
					{
						for(_ng1=0; _ng1<10; _ng1++)
						{
							BVF[_ng*10 + _ng1].data[j][i] = gfit_params_mpi[0].data[i][j][opt_gauss_number-1][_ng][_ng1];
							BVF_e[_ng*10 + _ng1].data[j][i] = gfit_params_e_mpi[0].data[i][j][opt_gauss_number-1][_ng][_ng1];
						}
					}
				}
			}
		}

		// 3. filtering weired rms level and update s/n
		int bn = 0;
		//for(i=0; i<TRparam[0].nax1; i++)
		for(i=xlower; i<xupper; i++)
		{
			//for(j=0; j<TRparam[0].nax2; j++)
			for(j=ylower; j<yupper; j++)
			{
				if(!isinf(raw_cube_mpi[0].data[i][j][0]) && !isnan(raw_cube_mpi[0].data[i][j][TRparam[0].nax3-1]) && !isinf(HI_VF[0].data[j][i]) && !isnan(HI_VF[0].data[j][i]))
				{
					// extract _boxfilter
					bn = 0;
					for(mi=-(box_x-1)/2; mi<(box_x+1)/2; mi++)
					{
						for(mj=-(box_y-1)/2; mj<(box_y+1)/2; mj++)
						{
							if(i+mi < 0 || i+mi >= TRparam[0].nax1 || j+mj < 0 || j+mj >= TRparam[0].nax2) continue;
							if(!isinf(BVF[6].data[j+mj][i+mi]) && !isnan(BVF[6].data[j+mj][i+mi]) && BVF[6].data[j+mj][i+mi] != 0) // if no blank
							{
								_filterbox[bn] = BVF[6].data[j+mj][i+mi];
								bn++;
							}
						}
					}
					robust_mean_std(_filterbox, bn, &hist_mean_filterbox, &hist_std_filterbox);
					//printf("%d %d %e %e %e\n", i, j, BVF[6].data[j][i], hist_mean_filterbox, hist_std_filterbox);

					if(BVF[6].data[j][i] < hist_mean_filterbox-5*hist_std_filterbox || BVF[6].data[j][i] > hist_mean_filterbox+5*hist_std_filterbox)
					{
						// update the s/n of each gaussian comps. w.r.t the new RMS
						for(_ng=0; _ng<opt_gauss_number; _ng++)
						{
							BVF[_ng*10 + 7].data[j][i] *= BVF[6].data[j][i]/hist_mean_filterbox;
							BVF_e[_ng*10 + 7].data[j][i] *= BVF[6].data[j][i]/hist_mean_filterbox;
						}

						// s/n of the highest: rescaling w.r.t the new RMS level
						BVF[8].data[j][i] *= BVF[6].data[j][i]/hist_mean_filterbox;
						// s/n of the lowest
						BVF[9].data[j][i] *= BVF[6].data[j][i]/hist_mean_filterbox;
						// new RMS
						BVF[6].data[j][i] = hist_mean_filterbox;
					}
				}
			}
		}


		// make directories if not exist
		sprintf(dirname, "%s/%s", argv[1], argv[9]);
		if(stat(dirname, &st) != 0)
		{
		    mkdir(dirname, 0775);
		}

		sprintf(dirname, "%s/%s/output%d", argv[1], argv[9], output_index);
		if(stat(dirname, &st) != 0)
		{
		    mkdir(dirname, 0775);
		}

		sprintf(dirname, "%s/%s/output%d", argv[1], argv[9], output_index);

		for(_ng=0; _ng<TRparam[0].n_gauss; _ng++)
		{
		    for(i=0; i<10; i++)
		    {
			sprintf(dirname, "%s/%s/output%d/G%02dg%02d", argv[1], argv[9], output_index, TRparam[0].n_gauss, _ng+1);
			if(stat(dirname, &st) != 0)
			{
			    mkdir(dirname, 0775);
			}
		    }
		}


		// 1. bulk directory
		sprintf(dirname, "%s/%s/output%d/bulk", argv[1], argv[9], output_index);
    	if(stat(dirname, &st) != 0)
		{
        	mkdir(dirname, 0775);
		}
		// 2. psgfit directory
		sprintf(dirname, "%s/%s/output%d/ps_gfit", argv[1], argv[9], output_index);
    	if(stat(dirname, &st) != 0)
		{
        	mkdir(dirname, 0775);
		}
		// 3. sgfit directory
		sprintf(dirname, "%s/%s/output%d/single_gfit", argv[1], argv[9], output_index);
    	if(stat(dirname, &st) != 0)
		{
        	mkdir(dirname, 0775);
		}
		// 4. primary_gfit directory
		sprintf(dirname, "%s/%s/output%d/primary_gfit", argv[1], argv[9], output_index);
    	if(stat(dirname, &st) != 0)
		{
        	mkdir(dirname, 0775);
		}
		// 5. cold gfit directory
		sprintf(dirname, "%s/%s/output%d/cold_gfit", argv[1], argv[9], output_index);
    	if(stat(dirname, &st) != 0)
		{
        	mkdir(dirname, 0775);
		}
		// 6. warm gfit directory
		sprintf(dirname, "%s/%s/output%d/warm_gfit", argv[1], argv[9], output_index);
    	if(stat(dirname, &st) != 0)
		{
        	mkdir(dirname, 0775);
		}
		// 7. strong nonc. gfit directory
		sprintf(dirname, "%s/%s/output%d/snonc_gfit", argv[1], argv[9], output_index);
    	if(stat(dirname, &st) != 0)
		{
        	mkdir(dirname, 0775);
		}
		// 8. weak nonc. gfit directory
		sprintf(dirname, "%s/%s/output%d/wnonc_gfit", argv[1], argv[9], output_index);
    	if(stat(dirname, &st) != 0)
		{
        	mkdir(dirname, 0775);
		}
		// 9. HVC gfit directory
		sprintf(dirname, "%s/%s/output%d/hvc_gfit", argv[1], argv[9], output_index);
    	if(stat(dirname, &st) != 0)
		{
        	mkdir(dirname, 0775);
		}


		_ng = 0;
		// 1. bulk
		for(i=0; i<10; i++)
		{
			sprintf(&fname[0].fitsfile_bvf[_ng][i], "%s/%s/output%d/bulk/%s.bulk.%d.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis[10*0+i]);

			sprintf(&fname[0].fitsfile_bvf_e[_ng][i], "%s/%s/output%d/bulk/%s.bulk.%d.e.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf_e[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis_e[10*0+i]);
		}

		// 2. psgfit directory
		for(i=0; i<10; i++)
		{
			sprintf(&fname[0].fitsfile_bvf[_ng][i], "%s/%s/output%d/ps_gfit/%s.ps_gfit.%d.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis[10*1+i]);

			sprintf(&fname[0].fitsfile_bvf_e[_ng][i], "%s/%s/output%d/ps_gfit/%s.ps_gfit.%d.e.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf_e[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis_e[10*1+i]);
		}

		// 3. sgfit directory
		for(i=0; i<10; i++)
		{
			sprintf(&fname[0].fitsfile_bvf[_ng][i], "%s/%s/output%d/single_gfit/%s.single_gfit.%d.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis[10*2+i]);

			sprintf(&fname[0].fitsfile_bvf_e[_ng][i], "%s/%s/output%d/single_gfit/%s.single_gfit.%d.e.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf_e[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis_e[10*2+i]);
		}

		// 4. primary_gfit directory
		for(i=0; i<10; i++)
		{
			sprintf(&fname[0].fitsfile_bvf[_ng][i], "%s/%s/output%d/primary_gfit/%s.primary_gfit.%d.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis[10*3+i]);

			sprintf(&fname[0].fitsfile_bvf_e[_ng][i], "%s/%s/output%d/primary_gfit/%s.primary_gfit.%d.e.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf_e[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis_e[10*3+i]);
		}

		// 5. cold gfit directory
		for(i=0; i<10; i++)
		{
			sprintf(&fname[0].fitsfile_bvf[_ng][i], "%s/%s/output%d/cold_gfit/%s.cold_gfit.%d.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis[10*4+i]);

			sprintf(&fname[0].fitsfile_bvf_e[_ng][i], "%s/%s/output%d/cold_gfit/%s.cold_gfit.%d.e.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf_e[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis_e[10*4+i]);
		}

		// 6. warm gfit directory
		for(i=0; i<10; i++)
		{
			sprintf(&fname[0].fitsfile_bvf[_ng][i], "%s/%s/output%d/warm_gfit/%s.warm_gfit.%d.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis[10*5+i]);

			sprintf(&fname[0].fitsfile_bvf_e[_ng][i], "%s/%s/output%d/warm_gfit/%s.warm_gfit.%d.e.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf_e[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis_e[10*5+i]);
		}

		// 7. strong nonc. gfit directory
		for(i=0; i<10; i++)
		{
			sprintf(&fname[0].fitsfile_bvf[_ng][i], "%s/%s/output%d/snonc_gfit/%s.snonc_gfit.%d.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis[10*6+i]);

			sprintf(&fname[0].fitsfile_bvf_e[_ng][i], "%s/%s/output%d/snonc_gfit/%s.snonc_gfit.%d.e.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf_e[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis_e[10*6+i]);
		}

		// 8. weak nonc. gfit directory
		for(i=0; i<10; i++)
		{
			sprintf(&fname[0].fitsfile_bvf[_ng][i], "%s/%s/output%d/wnonc_gfit/%s.wnonc_gfit.%d.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis[10*7+i]);

			sprintf(&fname[0].fitsfile_bvf_e[_ng][i], "%s/%s/output%d/wnonc_gfit/%s.wnonc_gfit.%d.e.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf_e[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis_e[10*7+i]);
		}

		// 9. HVC gfit directory
		for(i=0; i<10; i++)
		{
			sprintf(&fname[0].fitsfile_bvf[_ng][i], "%s/%s/output%d/hvc_gfit/%s.hvc_gfit.%d.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis[10*8+i]);

			sprintf(&fname[0].fitsfile_bvf_e[_ng][i], "%s/%s/output%d/hvc_gfit/%s.hvc_gfit.%d.e.fits", argv[1], argv[9], output_index, argv[2], i);
			save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf_e[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_analysis_e[10*8+i]);
		}


		// save all the slices
		for(_ng=0; _ng<TRparam[0].n_gauss; _ng++)
		{
			for(i=0; i<10; i++)
			{
				sprintf(&fname[0].fitsfile_bvf[_ng][i], "%s/%s/output%d/G%02dg%02d/%s.bvf.g%d.%d.fits", argv[1], argv[9], output_index, TRparam[0].n_gauss, _ng+1, argv[2], _ng, i);
				save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF[_ng*10 + i]);

				sprintf(&fname[0].fitsfile_bvf_e[_ng][i], "%s/%s/output%d/G%02dg%02d/%s.bvf.g%d.%d.e.fits", argv[1], argv[9], output_index, TRparam[0].n_gauss, _ng+1, argv[2], _ng, i);
				save_2dmapsfits(fptr_2d_refvf, fname[0].fitsfile_trfit_model, &fname[0].fitsfile_bvf_e[_ng][i], TRparam[0].nax1, TRparam[0].nax2, &BVF_e[_ng*10 + i]);
			}
		}

/*
i0 = xlower;
j0 = ylower;
//printf("\n\n%d %d sgfit: bg:%f a:%f s:%f x:%f %d sn:%f\n", i0, j0, gfit_params_mpi[0].data[i0][j0][0][0][0], gfit_params_mpi[0].data[i0][j0][0][0][1], gfit_params_mpi[0].data[i0][j0][0][0][2]/1000, gfit_params_mpi[0].data[i0][j0][0][0][3]/1000, (int)gfit_params_mpi[0].data[i0][j0][0][0][5], gfit_params_mpi[0].data[i0][j0][0][0][7]);
printf("\n\n%d %d sgfit: bg:%f a:%f s:%f x:%f %d sn:%f\n", i0, j0, gfit_params_mpi[0].data[i0][j0][0][0][0], gfit_params_mpi[0].data[i0][j0][0][0][1], gfit_params_mpi[0].data[i0][j0][0][0][2]/1000, gfit_params_mpi[0].data[i0][j0][0][0][3]/1000, (int)gfit_params_mpi[0].data[i0][j0][0][0][5], gfit_params_mpi[0].data[i0][j0][0][0][7]);
printf("\n\n%d %d sgfit: bge:%f ae:%f se:%f xe:%f %d sn:%f rms:%f\n", i0, j0, gfit_params_e_mpi[0].data[i0][j0][0][0][0], gfit_params_e_mpi[0].data[i0][j0][0][0][1], gfit_params_e_mpi[0].data[i0][j0][0][0][2]/1000, gfit_params_e_mpi[0].data[i0][j0][0][0][3]/1000, (int)gfit_params_mpi[0].data[i0][j0][0][0][5], gfit_params_mpi[0].data[i0][j0][0][0][7], gfit_params_mpi[0].data[i0][j0][0][0][6]);

//printf("\n\n%d %d 1: bg:%f a:%f s:%f x:%f %d sn:%f\n", i0, j0, gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][0], gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][1]/(sqrt(2.*M_PI)*gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][2]), gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][2]/1000, gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][3]/1000, (int)gfit_params_mpi[0].data[i0][j0][0][0][5], gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][7]);
printf("\n\n%d %d 1: bg:%f a:%f s:%f x:%f %d sn:%f\n", i0, j0, gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][0], gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][1], gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][2]/1000, gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][3]/1000, (int)gfit_params_mpi[0].data[i0][j0][0][0][5], gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][7]);
printf("\n\n%d %d 1: bge:%f ae:%f se:%f xe:%f %d sn:%f\n", i0, j0, gfit_params_e_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][0], gfit_params_e_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][1], gfit_params_e_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][2]/1000, gfit_params_e_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][3]/1000, (int)gfit_params_mpi[0].data[i0][j0][0][0][5], gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][7]);

//printf("\n\n%d %d 1: bg-e:%f a-e:%f s-e:%f x-e:%f %d sn:%f\n", i0, j0, gfit_params_e_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][0], gfit_params_e_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][1]/(sqrt(2.*M_PI)*gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][2]), gfit_params_e_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][2], gfit_params_e_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][0][3], (int)gfit_params_mpi[0].data[i0][j0][0][0][5], gfit_params_mpi[0].data[i0][j0][0][0][7]);

printf("%d %d 2: bg:%f a:%f s:%f x:%f %d\n", i0, j0, gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][1][0], gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][1][1], gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][1][2]/1000, gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][1][3]/1000, (int)gfit_params_mpi[0].data[i0][j0][0][0][5]);
printf("%d %d 2: bge:%f ae:%f se:%f xe:%f %d\n", i0, j0, gfit_params_e_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][1][0], gfit_params_e_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][1][1], gfit_params_e_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][1][2]/1000, gfit_params_e_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][1][3]/1000, (int)gfit_params_mpi[0].data[i0][j0][0][0][5]);
//printf("%d %d 2: bg:%f a:%f s:%f x:%f %d\n", i0, j0, gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][1][0], gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][1][1]/(sqrt(2.*M_PI)*gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][1][2]), gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][1][2]/1000, gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][1][3]/1000, (int)gfit_params_mpi[0].data[i0][j0][0][0][5]);

printf("%d %d 3: bg:%f a:%f s:%f x:%f %d\n", i0, j0, gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][2][0], gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][2][1], gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][2][2]/1000, gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][2][3]/1000, (int)gfit_params_mpi[0].data[i0][j0][0][0][5]);

printf("%d %d 4: bg:%f a:%f s:%f x:%f %d\n", i0, j0, gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][3][0], gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][3][1]/(sqrt(2.*M_PI)*gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][3][2]), gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][3][2]/1000, gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][3][3]/1000, (int)gfit_params_mpi[0].data[i0][j0][0][0][5]);

printf("%d %d 5: bg:%f a:%f s:%f x:%f %d\n", i0, j0, gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][4][0], gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][4][1]/(sqrt(2.*M_PI)*gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][4][2]), gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][4][2]/1000, gfit_params_mpi[0].data[i0][j0][(int)gfit_params_mpi[0].data[i0][j0][0][0][5]-1][4][3]/1000, (int)gfit_params_mpi[0].data[i0][j0][0][0][5]);

printf("low:%f high:%f std:%f\n", gfit_params_mpi[0].data[i0][j0][0][0][9], gfit_params_mpi[0].data[i0][j0][0][0][8], gfit_params_mpi[0].data[i0][j0][0][0][6]);
*/

		clock_t stop_time = clock();
		double elapsed = (double)(stop_time - start_time) / CLOCKS_PER_SEC;
		printf("THE EXECUTION TIME ELAPSED: %f SECS\n\n", elapsed);
	}

	free(_filterbox);
	free(vel_sort[0]);
	free(vel_sort);
	free(amp_sort[0]);
	free(amp_sort);
	free(disp_sort[0]);
	free(disp_sort);
	free(vel_gauss_n_sort[0]);
	free(vel_gauss_n_sort);
	free(HI_VF[0].data);
	free(HI_VF_boxFiltered[0].data);
	free(HI_VF_boxFiltered_sigma[0].data);

    //MPI_Finalized(&flag);
	MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

// --- End of line

