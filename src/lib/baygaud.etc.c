
#include "baygaud.etc.h"

// 2DBAT user defined functions
// ETC

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void read_user_input_baygaud(TR_ringParameters *TRparam, multinest_paramters *multinest_param, filename_2dbat *fname, int argc, char *argv[])
{
    int i=0;
    FILE *file_exist;

    // -----------------------------------------------------------------------
    // A. 2D maps ----------------------------
    // [1. WORKING DIRECTORY= /wdir]
    sprintf(TRparam[0].wdir, "%s", argv[1]);
    // [2. 3D data cube= cube.fits]
    sprintf(fname[0].raw_cube, "%s/%s", argv[1], argv[2]);


    // [3. 2D ref vf = model.vf.fits: use HI_VF global parameter for this]
    sprintf(fname[0].fitsfile_trfit_model, "%s/%s", argv[1], argv[3]);
    if((file_exist = fopen(fname[0].fitsfile_trfit_model, "r")) == NULL)
    {
		printf("%s doesn't exist... exit now\n", fname[0].fitsfile_trfit_model);
		exit(1);
    }

    TRparam[0].decim_x0 = 0;
    TRparam[0].decim_y0 = 0;
	TRparam[0]._nu_studenT = 3;

    // -----------------------------------------------------------------------
    //E. MULTINEST parameters : full run ----------------
    //multinest_param = (multinest_paramters *)malloc(sizeof (multinest_paramters) * 1);
    // [31. is= 0 or 1]
    multinest_param[0].is = 0;                 // important nested sampling?
    // [32. ceff= 0 or 1]
    multinest_param[0].ceff = 0;                   // run in constant efficiency mode?
    // [33. nlive= 50]
    multinest_param[0].nlive = 100;              // number of live points
    multinest_param[0].nlive_einasto_halofit = 100;              // number of live points
    // [34. efr= 0.8]
    multinest_param[0].efr = 0.8;                // set the required efficiency
    // [35. tol= 0.3]
    multinest_param[0].tol = 0.3;                // tol, defines the stopping criteria
    // [36. fb= 0]
    multinest_param[0].fb = 0;                 // need feedback on standard output?
    // [37. outfile= 0]
    multinest_param[0].outfile = 1;                // write output files?
    // [38. maxiter= 0]
    multinest_param[0].maxiter = 1000000;              // max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 

    //F. MULTINEST parameters : copy fullfit params ----------------
    multinest_param[0]._is = multinest_param[0].is;
    multinest_param[0]._ceff = multinest_param[0].ceff;
    multinest_param[0]._nlive_einasto_halofit = multinest_param[0].nlive_einasto_halofit;
    multinest_param[0]._efr = multinest_param[0].efr;
    multinest_param[0]._tol = multinest_param[0].tol;
    multinest_param[0]._fb = multinest_param[0].fb;
    multinest_param[0]._outfile = multinest_param[0].outfile;
    multinest_param[0]._maxiter = multinest_param[0].maxiter;

    // default multinest parameters
    multinest_param[0].mmodal = 1;                  // do mode separation?
    multinest_param[0].ndims = 999;                 // dimensionality (no. of free parameters)
    multinest_param[0].nPar = 999;                  // total no. of parameters including free & derived parameters
    multinest_param[0].nClsPar = 999;
    multinest_param[0].updInt = 200;                // after how many iterations feedback is required & the output files should be updated
                            // note: posterior files are updated & dumper routine is called after every updInt*10 iterations
    multinest_param[0].Ztol = -1E90;                // all the modes with logZ < Ztol are ignored
    multinest_param[0].maxModes = 100;              // expected max no. of modes (used only for memory allocation)
    for(i = 0; i < multinest_param[0].ndims; i++) multinest_param[0].pWrap[i] = 0;

    //strcpy(multinest_param[0].root, "%s/2dbat_output/multinest.");        // root for output files
    sprintf(multinest_param[0].root, "%s", argv[1]);        // root for output files
    multinest_param[0].seed = -1;       //          // random no. generator seed, if < 0 then take the seed from system clock
    multinest_param[0].resume = 0;                  // resume from a previous job?
    multinest_param[0].initMPI = 0;             // initialize MPI routines?, relevant only if compiling with MPI
                            // set it to F if you want your main program to handle MPI initialization
    multinest_param[0].logZero = -DBL_MAX;          // points with loglike < logZero will be ignored by MultiNest
                            // has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void usage_baygaud()
{
    printf("\n+ -------------------------------------------------------------------------------------- +\n");
    printf("+ Bayesian Automated GAussian Decomposer (BAGAD)                                         +\n");
    printf("+ by SE-HEON OH (KASI/ICRAR) + WALLABY KINEMATICS WORKING GROUP                          +\n");
    printf("+ -------------------------------------------------------------------------------------- +\n");
    printf("+                                                                                        +\n");
    printf("+ Development history                                                                    +\n");
    printf("+ : V.1.0.0 28/Feb/2018                                                                  +\n");
    printf("+                                                                                        +\n");
    printf("+ ! USAGE ------------------------------------------------------------------------------ +\n");
    printf("+                                                                                        +\n");
    printf("+: mpirun -np [0. N-threads= 8] ./baygaud                                                +\n");
	printf("+         [working-dir] [input-cube.fits] [input-ref-vf.fits]                            +\n");
	printf("+         [N-Gauss] [xlower] [ylower] [xupper] [yupper] [output_index(e.g., 0 ~ xxx)]    +\n");
    printf("+                                                                                        +\n");
    printf("+ A. INPUT data ------------------------------------------------------------------------ +\n");
    printf("+  [1. WORKING DIRECTORY= /home/baygaud/wdir]                                            +\n");
	printf("+    : where the [input-cube.fits] and [input-ref-vf.fits] locates                       +\n");
    printf("+  [2. INPUT DATA CUBE= input-cube.fits]                                                 +\n");
    printf("+  [3. INPUT 2D REFERENCE VELOITY FIELD= input-ref-vf.fits in (km/s)]                    +\n");
    printf("+  ..................................................................................... +\n");
    printf("+ B. The maximum number of Gaussian components to be decomposed ------------------------ +\n");
    printf("+  [4. N-Gauss= 4]                                                                       +\n");
    printf("+  ..................................................................................... +\n");
    printf("+ C. AREA to extract ------------------------------------------------------------------- +\n");
    printf("+  [5. xlower= 0 in pixel]                                                               +\n");
    printf("+  [6. ylower= 0 in pixel]                                                               +\n");
    printf("+  [7. xupper= nax1 in pixel]                                                            +\n");
    printf("+  [8. yupper= nax2 in pixel]                                                            +\n");
    printf("+  ..................................................................................... +\n");
    printf("+ D. OUTPUT directory index ------------------------------------------------------------ +\n");
    printf("+  [9. output_index= 0 or whatever any number]                                           +\n");
    printf("+  ------------------------------------------------------------------------------------- +\n");
    printf("+                                                                                        +\n");
    printf("+ ! EXAMPLE ---------------------------------------------------------------------------- +\n");
    printf("+                                                                                        +\n");
    printf("+ mpirun -np 4 ./baygaud \\                                                               +\n");
    printf("+ /home/seheon/baygaud.v1.0.0/baygaud/ngc2403 \\                                          +\n");
    printf("+ input-cube.fits \\                                                                      +\n");
    printf("+ input-ref-vf.fits \\                                                                    +\n");
    printf("+ 4 \\                                                                                    +\n");
    printf("+ 100 100 140 140 \\                                                                      +\n");
    printf("+ 0                                                                                      +\n");
    printf("+                                                                                        +\n");
    printf("+ -------------------------------------------------------------------------------------- +\n\n");
}



// --- End of line

