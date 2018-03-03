#include "baygaud.gfit.h"

// 2DBAT user defined functions
// Gaussian fit related


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Gfit_multinest(multinest_paramters *multinest_param, TR_ringParameters *TRparam)
{
    /* set the MultiNest sampling parameters */
    char root[500];     // root for output files
    int i=0;

    int is, mmodal, ceff, nlive, ndims, nPar, nClsPar, updInt, maxModes, seed, fb, resume, outfile, initMPI, maxiter;
    double efr, tol, Ztol, logZero;

    /* set the MultiNest sampling parameters */
    is = multinest_param[0].is;
    mmodal = multinest_param[0].mmodal;
    ceff = multinest_param[0].ceff;
    nlive = multinest_param[0].nlive;
    efr = multinest_param[0].efr;
    tol = multinest_param[0].tol;
    updInt = multinest_param[0].updInt;

    Ztol = multinest_param[0].Ztol;
    maxModes = multinest_param[0].maxModes;
    strcpy(root, "gfit.");

    seed = multinest_param[0].seed;
    fb = multinest_param[0].fb;
    resume = multinest_param[0].resume;
    //outfile = multinest_param[0].outfile;
    outfile = 0; // for avoiding name output name confliction between different ranks
    initMPI = multinest_param[0].initMPI;
    logZero = multinest_param[0].logZero;
    maxiter = multinest_param[0].maxiter;

    ndims = 3*TRparam[0].n_gauss+1;
    nPar = 3*TRparam[0].n_gauss+1;
    nClsPar = 3*TRparam[0].n_gauss+1;

    int pWrap[ndims];
    for(i=0; i<ndims; i++)
    {
        pWrap[i] = multinest_param[0].pWrap[i];
    }

    /* Calling multinest */
    run_mpioff(is, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, loglikelihood_gfit, dumper_Gfits, TRparam);
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Gfit_multinest_student(multinest_paramters *multinest_param, TR_ringParameters *TRparam)
{
    /* set the MultiNest sampling parameters */
    char root[500];     // root for output files
    int i=0;

    int is, mmodal, ceff, nlive, ndims, nPar, nClsPar, updInt, maxModes, seed, fb, resume, outfile, initMPI, maxiter;
    double efr, tol, Ztol, logZero;

    /* set the MultiNest sampling parameters */
    is = multinest_param[0].is;
    mmodal = multinest_param[0].mmodal;
    ceff = multinest_param[0].ceff;
    nlive = multinest_param[0].nlive;
    efr = multinest_param[0].efr;
    tol = multinest_param[0].tol;
    updInt = multinest_param[0].updInt;

    Ztol = multinest_param[0].Ztol;
    maxModes = multinest_param[0].maxModes;
    strcpy(root, "gfit.");

    seed = multinest_param[0].seed;
    fb = multinest_param[0].fb;
    resume = multinest_param[0].resume;
    //outfile = multinest_param[0].outfile;
    outfile = 0; // for avoiding name output name confliction between different ranks
    initMPI = multinest_param[0].initMPI;
    logZero = multinest_param[0].logZero;
    maxiter = multinest_param[0].maxiter;

    ndims = 3*TRparam[0].n_gauss+1+1; // +1 for e_sigma_fitted...
    nPar = 3*TRparam[0].n_gauss+1+1; // +1 for e_sigma_fitted...
    nClsPar = 3*TRparam[0].n_gauss+1+1; // +1 for e_sigma_fitted...

    int pWrap[ndims];
    for(i=0; i<ndims; i++)
    {
        pWrap[i] = multinest_param[0].pWrap[i];
    }

    /* Calling multinest */
    run_mpioff(is, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, loglikelihood_gfit_student, dumper_Gfits_student, TRparam);
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void loglikelihood_gfit(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam)
{
    int i=0, j=0;
    int n_gauss;

    double chi2 = 0.0;
    double logsum_errors = 0.;
    double GLLhood0 = 0.0;
    double slhood = 0.0;
    double npoints = 0.0;
    double gauss_model=0;
	double vel=0, vel_norm=0, vel_min=0, vel_max=0, vel_temp;
	double _r, _mu, _sigma;

	n_gauss = TRparam[0].n_gauss;

    /* set uniform priors x1 ~ x2*/
    /* convert unit Cube to actual parameter values */
    // 1. priors for coefficients: x1 ~ x2


    for(i=0; i<n_gauss; i++)
    {
        if(i==0)
        {
            Cube[i] = TRparam[0].g01 + Cube[i]*(TRparam[0].g02-TRparam[0].g01); // coefficient g0
            Cube[1+3*i] = TRparam[0].gA1 + Cube[1+3*i]*(TRparam[0].gA2-TRparam[0].gA1); // coefficients A
            Cube[2+3*i] = TRparam[0].gS1 + Cube[2+3*i]*(TRparam[0].gS2-TRparam[0].gS1); // coefficients Sigma
            Cube[3+3*i] = TRparam[0].gX1 + Cube[3+3*i]*(TRparam[0].gX2-TRparam[0].gX1); // coefficients X
        }
        else
        {
            Cube[1+3*i] = TRparam[0].gA1 + Cube[1+3*i]*(TRparam[0].gA2-TRparam[0].gA1); // coefficients A
            Cube[2+3*i] = TRparam[0].gS1 + Cube[2+3*i]*(TRparam[0].gS2-TRparam[0].gS1); // coefficients Sigma
            Cube[3+3*i] = TRparam[0].gX1 + Cube[3+3*i]*(TRparam[0].gX2-TRparam[0].gX1); // coefficients X
        }
    }

/*
    for(i=0; i<n_gauss; i++)
    {
        if(i==0)
        {
			// bg
        	_r = Cube[i];
        	_mu = 0.5;
        	_sigma = 0.1;
            Cube[i] = fabs(gaussian_prior(_r, _mu, _sigma));

			// Amp
        	_r = Cube[1+3*i];
        	_mu = 0.1;
        	_sigma = 0.03;
            Cube[1+3*i] = fabs(gaussian_prior(_r, _mu, _sigma));

			// sigma
        	_r = Cube[2+3*i];
        	_mu = 0.2;
        	_sigma = 0.1; 
            Cube[2+3*i] = gaussian_prior(_r, _mu, _sigma);

			// X
        	_r = Cube[3+3*i];
        	_mu = 0.5;
        	_sigma = 0.17;
            Cube[3+3*i] = gaussian_prior(_r, _mu, _sigma);
        }
        else
        {
			// Amp
        	_r = Cube[1+3*i];
        	_mu = 0.1;
        	_sigma = 0.03;
            Cube[1+3*i] = fabs(gaussian_prior(_r, _mu, _sigma));

			// sigma
        	_r = Cube[2+3*i];
        	_mu = 0.2;
        	_sigma = 0.1; 
            Cube[2+3*i] = gaussian_prior(_r, _mu, _sigma);

			// X
        	_r = Cube[3+3*i];
        	_mu = 0.5;
        	_sigma = 0.17;
            Cube[3+3*i] = gaussian_prior(_r, _mu, _sigma);
        }
    }
*/


    // 8. Total number of n_hist 
	npoints = (double)(TRparam[0].effective_channel_end-TRparam[0].effective_channel_start);
	GLLhood0 = -(npoints/2.0)*log(2.0*M_PI); // for individual error


	vel_min = (0 - TRparam[0].velocity_offset - TRparam[0].nax3)*TRparam[0].Velocity_Seperation + TRparam[0].Reference_Velocity;
	vel_max = (3*TRparam[0].nax3-1 - TRparam[0].velocity_offset - TRparam[0].nax3)*TRparam[0].Velocity_Seperation + TRparam[0].Reference_Velocity;

	if(vel_min > vel_max)
	{
		vel_temp = vel_min;
		vel_min = vel_max;
		vel_max = vel_temp;
	}

    chi2 = 0.;
    logsum_errors = 0.;
	for(i=0; i<3*TRparam[0].nax3; i++)
    {
		vel = (i - TRparam[0].velocity_offset - TRparam[0].nax3)*TRparam[0].Velocity_Seperation + TRparam[0].Reference_Velocity;
		vel_norm = (vel-vel_min)/(vel_max-vel_min); // normalise velocity to make it easier to fit

		//TRparam[0].velocity_profile_e[i] = 0.1;
		gauss_model = gauss_function(Cube, n_gauss, vel_norm); // calculate a model value for a given polynomial function
		logsum_errors += log(TRparam[0].velocity_profile_e[i]);

//printf("%d %f %f %f %f\n", i, vel_norm, TRparam[0].velocity_profile_norm[i], gauss_model, TRparam[0].velocity_profile_e[i]);

		chi2 += powf(((TRparam[0].velocity_profile_norm[i] - gauss_model)/powf(TRparam[0].velocity_profile_e[i], 1)), 2);
	}
	slhood = GLLhood0 - logsum_errors - chi2/2.0;
	//slhood = GLLhood0 - chi2/2.0;
	*lnew = slhood;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void loglikelihood_gfit_student(double *Cube, int *ndim, int *npars, double *lnew, TR_ringParameters *TRparam)
{
    int i=0, j=0;
    int n_gauss;

    double chi2 = 0.0;
    double GLLhood0 = 0.0;
    double slhood = 0.0;
    double npoints = 0.0;
    double gauss_model=0;
	double vel=0, vel_norm=0, vel_min=0, vel_max=0, vel_temp;
	double _r, _mu, _sigma;

    // for student-t distribution parameters
	double logsum_sigma2=0, logsum_student_chi=0, _nu, log_Likelihood_studenT;
	double NobsPoints_available;
	double e_sigma;
	double logsum_errors = 0., sum_errors=0;


	n_gauss = TRparam[0].n_gauss;

    /* set uniform priors x1 ~ x2*/
    /* convert unit Cube to actual parameter values */
    // 1. priors for coefficients: x1 ~ x2


    for(i=0; i<n_gauss; i++)
    {
        if(i==0)
        {
            Cube[i] = TRparam[0].g01 + Cube[i]*(TRparam[0].g02-TRparam[0].g01); // coefficient g0
            Cube[1+3*i] = TRparam[0].gA1 + Cube[1+3*i]*(TRparam[0].gA2-TRparam[0].gA1); // coefficients A
            Cube[2+3*i] = TRparam[0].gS1 + Cube[2+3*i]*(TRparam[0].gS2-TRparam[0].gS1); // coefficients Sigma
            Cube[3+3*i] = TRparam[0].gX1 + Cube[3+3*i]*(TRparam[0].gX2-TRparam[0].gX1); // coefficients X
        }
        else
        {
            Cube[1+3*i] = TRparam[0].gA1 + Cube[1+3*i]*(TRparam[0].gA2-TRparam[0].gA1); // coefficients A
            Cube[2+3*i] = TRparam[0].gS1 + Cube[2+3*i]*(TRparam[0].gS2-TRparam[0].gS1); // coefficients Sigma
            Cube[3+3*i] = TRparam[0].gX1 + Cube[3+3*i]*(TRparam[0].gX2-TRparam[0].gX1); // coefficients X
        }
    }

	//Cube[3+3*i] = TRparam[0].gX1 + Cube[3+3*i]*(TRparam[0].gX2-TRparam[0].gX1); // coefficients X
	Cube[3*n_gauss+1] = 0.0 + Cube[3*n_gauss+1]*(0.03-0.0); // n_free_params = 3*n_gauss + 1 (0.1 ~ 0.3) : given the flux (0 - 1.0)
	e_sigma = Cube[3*n_gauss + 1];

	//printf("%f\n", e_sigma);

////////////////////////////////////////////////////////////////////


	vel_min = (0 - TRparam[0].velocity_offset)*TRparam[0].Velocity_Seperation + TRparam[0].Reference_Velocity;
	vel_max = (1*TRparam[0].nax3-1 - TRparam[0].velocity_offset)*TRparam[0].Velocity_Seperation + TRparam[0].Reference_Velocity;
	if(vel_min > vel_max)
	{
		vel_temp = vel_min;
		vel_min = vel_max;
		vel_max = vel_temp;
	}

    // 8. calculate log likelihood value
    chi2 = 0.;
    logsum_errors = 0.;
    sum_errors = 0.;
	//NobsPoints_available = (double)(TRparam[0].effective_channel_end-TRparam[0].effective_channel_start);

    // _nu of student-T distribution : _nu = 30 for normal : _nu = 1 recommended for best removing the outliers
    _nu = TRparam[0]._nu_studenT;
    _nu = 3;
    logsum_student_chi = 0;
    logsum_sigma2 = 0;

	NobsPoints_available = 0;
	for(i=0; i<1*TRparam[0].nax3; i++)
    {
		vel = (i - TRparam[0].velocity_offset)*TRparam[0].Velocity_Seperation + TRparam[0].Reference_Velocity;
		vel_norm = (vel-vel_min)/(vel_max-vel_min); // normalise velocity to make it easier to fit

		gauss_model = gauss_function(Cube, n_gauss, vel_norm); // calculate a model value for a given polynomial function


		NobsPoints_available += 1;
        logsum_errors += 0.5*log(e_sigma*e_sigma);
        logsum_student_chi += ((1+_nu)/2.0)*log(1.0+pow((TRparam[0].velocity_profile_norm[i] - gauss_model)/e_sigma, 2)/(_nu-2));
    }

//log_Likelihood_studenT = NobsPoints_available*(log(gsl_sf_gamma((_nu+1)/2.0)) - log(gsl_sf_gamma(_nu/2.0))) - logsum_errors - logsum_student_chi;
    log_Likelihood_studenT = NobsPoints_available*(log(gsl_sf_gamma((_nu+1)/2.0)) - log(sqrt(M_PI*(_nu-2))*gsl_sf_gamma(_nu/2.0))) - logsum_errors - logsum_student_chi;
    *lnew = log_Likelihood_studenT;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void dumper_Gfits(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, TR_ringParameters *TRparam)
{
    // convert the 2D Fortran arrays to C arrays
    int i, j;
	int _i_, _j_;

    // the posterior distribution
    // postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
    double postdist[*nSamples][*nPar + 2];
    for( i = 0; i < *nPar + 2; i++ )
    {
        for( j = 0; j < *nSamples; j++ )
        {
            postdist[j][i] = posterior[0][i * (*nSamples) + j];
        }

    }

    // last set of live points
    // pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
    double pLivePts[*nlive][*nPar + 1];
    for( i = 0; i < *nPar + 1; i++ )
    {
        for( j = 0; j < *nlive; j++ )
        {
            pLivePts[j][i] = physLive[0][i * (*nlive) + j];
        }
    }

	// best fits
	//_i_ = TRparam[0]._i_;
	//_j_ = TRparam[0]._j_;

	TRparam[0].g_param[0] = paramConstr[0][*nPar*2+0]; // background first
	TRparam[0].g_param[TRparam[0].n_gauss*3+1] = paramConstr[0][*nPar*1+0]; // background std first

	for(j=1; j<3*TRparam[0].n_gauss+1; j++)
	{
        TRparam[0].g_param[j] = paramConstr[0][*nPar*2+j];
	}
	for(j=1; j<3*TRparam[0].n_gauss+1; j++)
	{
        TRparam[0].g_param[j+3*TRparam[0].n_gauss+1] = paramConstr[0][*nPar*1+j];
	}
			
    TRparam[0].maxLogLikeF = *maxLogLike;
    //TRparam[0].logZF = *logZ;
    //TRparam[0].logZerrF = *logZerr;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void dumper_Gfits_student(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *INSlogZ, double *logZerr, TR_ringParameters *TRparam)
{
    // convert the 2D Fortran arrays to C arrays
    int i, j;
	int _i_, _j_;
	int n_freeParams;

    // the posterior distribution
    // postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
    double postdist[*nSamples][*nPar + 2];
    for( i = 0; i < *nPar + 2; i++ )
    {
        for( j = 0; j < *nSamples; j++ )
        {
            postdist[j][i] = posterior[0][i * (*nSamples) + j];
        }

    }

    // last set of live points
    // pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
    double pLivePts[*nlive][*nPar + 1];
    for( i = 0; i < *nPar + 1; i++ )
    {
        for( j = 0; j < *nlive; j++ )
        {
            pLivePts[j][i] = physLive[0][i * (*nlive) + j];
        }
    }

	// best fits
	//_i_ = TRparam[0]._i_;
	//_j_ = TRparam[0]._j_;

	TRparam[0].g_param[0] = paramConstr[0][*nPar*2+0]; // background first
	TRparam[0].g_param[TRparam[0].n_gauss*3+1] = paramConstr[0][*nPar*1+0]; // background std first

	n_freeParams = 1;
	for(j=1; j<3*TRparam[0].n_gauss+1; j++)
	{
        TRparam[0].g_param[j] = paramConstr[0][*nPar*2+j];
		n_freeParams++;
	}
	for(j=1; j<3*TRparam[0].n_gauss+1; j++)
	{
        TRparam[0].g_param[j+3*TRparam[0].n_gauss+1] = paramConstr[0][*nPar*1+j];
	}

	TRparam[0].e_sigma_fitted = paramConstr[0][*nPar*2+n_freeParams]; // the fitted rms (std) of a current profile 
	//TRparam[0].e_sigma_fitted_error = paramConstr[0][*nPar*2+n_freeParams]; // not yet implemented...
    TRparam[0].maxLogLikeF = *maxLogLike;
    //TRparam[0].logZF = *logZ;
    //TRparam[0].logZerrF = *logZerr;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gauss_function(double *Cube, int n_gauss, double x)
{
    int i=0;
    double gauss_model = 0.;

    // G0 = Cube[0] : background
    // A0 = Cube[1+3*n_gauss], sigma = Cube[2+3*n_gauss], x0 = Cube[3+3*n_gauss];

    gauss_model += Cube[0];
    for(i=0; i<n_gauss; i++)
    {
        gauss_model += (Cube[1+3*i]/(sqrt(2.0*M_PI)*Cube[2+3*i])) * exp(-0.5*pow((x-Cube[3+3*i])/Cube[2+3*i], 2));
    }

    return gauss_model;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double gauss_function_at_a_pixel_gfit_params_mpi(float *****gfit_params, int n_gauss, int _i_, int _j_, double x)
{
    int i=0, opt_n_gauss;
    double gauss_model = 0.;

    // G0 = g_param[_i_][_j_][0] : background
    // A0 = g_param[_i_][_j_][1+3*n_gauss], sigma = g_param[_i_][_j_][2+3*n_gauss], x0 = g_param[_i_][_j_][3+3*n_gauss];

	//opt_n_gauss = (int)gfit_params[_i_][_j_][0][0][5];
    //gauss_model += gfit_params[_i_][_j_][opt_n_gauss-1][0][0];

	opt_n_gauss = (int)gfit_params[_i_][_j_][0][0][5];
    gauss_model += gfit_params[_i_][_j_][n_gauss-1][0][0];

    for(i=0; i<n_gauss; i++)
    {
        //gauss_model += (g_param[1+3*i]/(sqrt(2.0*M_PI)*g_param[2+3*i])) * exp(-0.5*pow((x-g_param[3+3*i])/g_param[2+3*i], 2));
        gauss_model += (gfit_params[_i_][_j_][n_gauss-1][i][1]/(sqrt(2.0*M_PI)* gfit_params[_i_][_j_][n_gauss-1][i][2] )) * exp(-0.5*pow((x- gfit_params[_i_][_j_][n_gauss-1][i][3] )/gfit_params[_i_][_j_][n_gauss-1][i][2], 2));
    }

    return gauss_model;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double each_gauss_function_at_a_pixel_gfit_params_mpi(float *****gfit_params, int n_gauss, int nth_gauss, int _i_, int _j_, double x)
{
    int i=0;
    double gauss_model = 0.;

    // G0 = g_param[_i_][_j_][0] : background
    // A0 = g_param[_i_][_j_][1+3*n_gauss], sigma = g_param[_i_][_j_][2+3*n_gauss], x0 = g_param[_i_][_j_][3+3*n_gauss];

	//opt_n_gauss = (int)gfit_params[_i_][_j_][0][0][5];
    gauss_model += gfit_params[_i_][_j_][n_gauss-1][0][0];

    //gauss_model = (gfit_params[_i_][_j_][n_gauss-1][nth_gauss][1]/(sqrt(2.0*M_PI)* gfit_params[_i_][_j_][n_gauss-1][nth_gauss][2] )) * exp(-0.5*pow((x- gfit_params[_i_][_j_][n_gauss-1][nth_gauss][3] )/gfit_params[_i_][_j_][n_gauss-1][nth_gauss][2], 2));
    gauss_model = (gfit_params[_i_][_j_][n_gauss-1][nth_gauss][1]/(sqrt(2.0*M_PI)* gfit_params[_i_][_j_][n_gauss-1][nth_gauss][2] ));

    return gauss_model;
}



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void ext_bulk_motions(int xlower, int xupper, int ylower, int yupper, multinest_paramters *multinest_param, TR_ringParameters *TRparam, int sorting_order, float ***HI_cube_mpi, int rank)
{
	int i=0, j=0, k=0, n=0, i_max, j_max, b0=0;
	int i_std=0;
	int status=0, anynul=0, group=0;
	int N_perfect_single_gaussian=0;
	int calc_initial_sigma=0;
	int nax1, nax2, nax3;
	int n_gauss_upto=0, n_gauss_upto_t=0;
	int _n0=0, _n1=0, _n2=0, _n3=0, _p=0;
	int _i_, _j_;

	nax1 = TRparam[0].nax1;
	nax2 = TRparam[0].nax2;
	nax3 = TRparam[0].nax3;

	double total_N;
	double sGfit_flux, dGfit_flux;
	double gauss_model, gauss_model0, gauss_model1, gauss_model2;
	double unit_channel_width=0, N_vel_data;
	double single_gauss_model;

	/* set ring parameters */
	double maxLogLike_sGfit, n_freeParams_sGfit, maxLogLike_dGfit, n_freeParams_dGfit;
	double vel, vel_norm, vel_min, vel_max, amp_min, amp_max, vel_temp;
	double gfit_BIC=0;
	double BIC_sGfit=0., BIC_dGfit=0., ratio_sGfit_dGfit_BIC=0;
	//double vel_profile_res[nax3];
	//double *velocity_profile_temp = (double*)malloc(sizeof(double) * nax3); // 1D array: Nrings x 1
	double velocity_profile_temp[3*nax3];
	double vel_profile_res_RBM=0., vel_profile_res_SIGMA=0.;
	double perfect_sGfit_amp=0., perfect_sGfit_sigma=0.;
	double bg_left_edge_rbm=0, bg_left_edge_sigma=0;
	double gauss1_flux, gauss2_flux;

	// 1. the represenative profile
	double max_profile_g01, max_profile_g02;
	double max_profile_gA1, max_profile_gA2;
	double max_profile_gS1, max_profile_gS2;
	double max_profile_gX1, max_profile_gX2;

	double peak_flux, peak_flux_t;

	int std_n;
	double std, sn_profile, _Xo;
	double res_mean, res_std;
	double n_total_profiles=0, n_processed_profiles=0;

	double Xo_lower_std, Xo_upper_std;

	double a_norm, s_norm, x_norm, g0_norm;
	double ae_norm, se_norm, xe_norm, g0e_norm;

	/* initialize dynamic 3D array */
	//double amp_profile[NAX2][NAX1][2];
	//double amp_profile_mpi[2][NAX1][NAX2];

	//double amp_profile[nax2][nax1][2];
	//double **amp_profile = (double**)malloc(sizeof(double *) * (int)((TRparam[0].xupper-TRparam[0].xlower)/TRparam[0].decimX)*(int)((TRparam[0].yupper-TRparam[0].ylower)/TRparam[0].decimY)); // 3D array: N_total x 3 
	double amp_profile_min, amp_profile_max;
	total_N = (int)((TRparam[0].xupper-TRparam[0].xlower)/TRparam[0].decimX)*(int)((TRparam[0].yupper-TRparam[0].ylower)/TRparam[0].decimY); //total number of profiles in the given region

    /* initialize dynamic 2D array */
    double **gfit_bic_temp = (double**)malloc(sizeof(double *) * TRparam[0].n_gauss); // 2D array: TRparam[0].n_gauss x 2
    for(i=0; i<TRparam[0].n_gauss; i++)
    {
       gfit_bic_temp[i] = (double*)malloc(sizeof(double) * 2);
    }

    /* initialize dynamic 2D array */
    double **gfit_A_temp = (double**)malloc(sizeof(double *) * TRparam[0].n_gauss); // 2D array: TRparam[0].n_gauss x 2
    for(i=0; i<TRparam[0].n_gauss; i++)
    {
       gfit_A_temp[i] = (double*)malloc(sizeof(double) * 2);
    }

	vel_min = (0 - TRparam[0].velocity_offset)*TRparam[0].Velocity_Seperation + TRparam[0].Reference_Velocity;
	vel_max = (1*nax3-1 - TRparam[0].velocity_offset)*TRparam[0].Velocity_Seperation + TRparam[0].Reference_Velocity;
	if(vel_min > vel_max)
	{
		vel_temp = vel_min;
		vel_min = vel_max;
		vel_max = vel_temp;
	}

	TRparam[0].unit_channel_width = 1.0/(double)(nax3*1);
	// BULK LIMIT 
	//bulk_limit = 0.2; // in normalised unit

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// I. Extract a velocity field which only has "perfect" single Gaussian velocity profile.
	// : A strict BIC criteria based small errors in velocity profile is used,
	// : i.e., TRparam[0].velocity_profile_e[i] = TRparam[0].vel_profile_res_SIGMA*0.2;

	n_total_profiles = (xupper-xlower)*(yupper-ylower);
	n_gauss_upto = TRparam[0].n_gauss;
	n_gauss_upto_t = TRparam[0].n_gauss;

	for(i=xlower; i<xupper; i++) // in the given region : xlower compared to TRparam[0].xlower
	{
		for(j=ylower; j<yupper; j++)
		{
			n_processed_profiles += 1;
			TRparam[0]._i_ = i; // pixel xcoord i
			TRparam[0]._j_ = j; // pixel ycoord j
			//if(HI_cube_mpi[i][j][(int)(nax3/2)] + 1 > HI_cube_mpi[i][j][(int)(nax3/2)])
			if(!isnan(HI_cube_mpi[i][j][0]) && !isinf(HI_cube_mpi[i][j][nax3-1]))
			//if(HI_VF[j][i] + 1 > HI_VF[j][i]) // if the current pixel is not blank in mom1
			{
				for(k=0; k<nax3; k++)
				{
					if(HI_cube_mpi[i][j][k]+1 > HI_cube_mpi[i][j][k]) // not blank	
					{
						velocity_profile_temp[k] = HI_cube_mpi[i][j][k]*1E5;
					}
					else
					{
						velocity_profile_temp[k] = 0;
					}
					//printf("%d %f\n", k, velocity_profile_temp[k]);
				}
				// sort the velocity profile to estimiate its maximum amplitude which is used for calculating an initial estimate for error
				qsort(&velocity_profile_temp[0], nax3, sizeof(double), cmp);
				for(k=0; k<nax3; k++)
				{
					velocity_profile_temp[k] = velocity_profile_temp[k]/1E5;
				}
				amp_profile_min = velocity_profile_temp[nax3-1]; // minimum amp
				amp_profile_max = velocity_profile_temp[0]; // maximum amp
				TRparam[0].vel_profile_res_SIGMA = velocity_profile_temp[0]*0.01; // set profile errors to 1% of the maximum amplitude 
				//free(velocity_profile_temp);

				// normalise profile
				for(k=0; k<1*nax3; k++)
				{
					if(HI_cube_mpi[i][j][k]+1 > HI_cube_mpi[i][j][k]) // not blank	
					{
						TRparam[0].velocity_profile[k] = HI_cube_mpi[i][j][k];
					}
					else
					{
						TRparam[0].velocity_profile[k] = 1E-5;
					}
					TRparam[0].velocity_profile_norm[k] = (TRparam[0].velocity_profile[k]-amp_profile_min)/(amp_profile_max-amp_profile_min);
				}

				// NOTE: adjust the priors (e.g., sigma values) and the tolerance level : 0.1 would be fine..
				//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				// 1. single Gaussian fit
				/* update uniform priors for gaussian function based on the above single Gaussian fit */
				// 1. background
				TRparam[0].g01 = 0.0;
				TRparam[0].g02 = 0.3;
				// 2. sigma
				//TRparam[0].gS1 = 0.01*fabs(TRparam[0].Velocity_Seperation / (vel_max - vel_min)); // channel resolution in unit
				TRparam[0].gS1 = 0.01;
				TRparam[0].gS2 = 0.2;
				// 3. amplitude
				//TRparam[0].gA1 = 0.01 * sqrt(2*M_PI)*TRparam[0].gS1; // this is calculated from sqrt(2pi)*gS1
				//TRparam[0].gA1 = 0.01*sqrt(2*M_PI)*TRparam[0].gS1;
				TRparam[0].gA1 = 0.0;
				//TRparam[0].gA2 = 1.0 * sqrt(2*M_PI)*TRparam[0].gS2; // this is calculated from sqrt(2pi)*gS2
				TRparam[0].gA2 = 0.1;

				// 4. centre
				TRparam[0].gX1 = 0.1;
				TRparam[0].gX2 = 0.9;

				// start multinest
				multinest_param[0].nlive = 80;
				multinest_param[0].maxiter = 100000;
				multinest_param[0].tol = 0.1;
				multinest_param[0].efr = 0.9;

				// set effective channels
				//TRparam[0].effective_channel_start = (int)(TRparam[0].gX1/TRparam[0].unit_channel_width);
				//TRparam[0].effective_channel_end = (int)(TRparam[0].gX2/TRparam[0].unit_channel_width);
				TRparam[0].effective_channel_start = 0;
				TRparam[0].effective_channel_end = 1*nax3;
				N_vel_data = (double)(TRparam[0].effective_channel_end - TRparam[0].effective_channel_start);

					
				// ++++++++++++++++++++++++++++++++++++++++++++++++++++
				// gfit_params_mpi description

				// gfit_params_mpi[i][j]  [gauss-#:0, 1, 2...][gauss-sequence-#:0, 1, 2...][param-#:0(bg), 1(A), 2(sig), 3(X), 4(BIC), 5(optimal gaussn-#] 6(std), 7(s/n of each gaussian components)], 8[s/n of the HIGEST FLUX] 9[s/n of the LOWEST gaussian component]: note the HIGHEST FLUX means not the fitteted Gaussian component with the highest flux but the total peak flux

				// gfit_params_mpi[i][j]  [gauss-#:0, 1, 2...][gauss-sequence-#:0, 1, 2...][param-#:0(bg), 1(A), 2(sig), 3(X), 4(BIC), 5(perfect single Gauss] 6(std), 7(s/n)], 8[optimal Gauss-#] 9[s/n same as in 7?]
				// gfit_params_e_mpi[i][j][gauss-#:0, 1, 2...][gauss-sequence-#:0, 1, 2...][param-#:0(bg_e), 1(A_e), 2(sig_e), 3(X_e), 4(none), 5(none), 6(none), 7(none)]
				// gfit_params_mpi[i][j][0][0][5] = optimal gauss-# based on BIC
				// note gauss-# starts from 0: e.g., if the optimal gauss-# found is 3, then gfit_params_mpi[i][j][3-1][0][0], ....
				// ++++++++++++++++++++++++++++++++++++++++++++++++++++

				

				for(i_std=0; i_std<1; i_std++)
				{
					//if(i_std==0)
					//{
						//for(k=0; k<1*nax3; k++) // start with initial error values...
						//{
						//	TRparam[0].velocity_profile_e[k] = 1E-2; // default uncertainties for single gauss fit
						//}
					//}

					// initialize
					for(k=0; k<n_gauss_upto_t; k++)
					{
					   gfit_bic_temp[k][0] = 0;
					   gfit_bic_temp[k][1] = 1E9; // put a large value
					}

					for(_n0=1; _n0<n_gauss_upto+1; _n0++)
					{
						TRparam[0].n_gauss = _n0;
						Gfit_multinest_student(multinest_param, TRparam);

						for(_n1=0; _n1<_n0; _n1++)
						{
							gfit_params_mpi[0].data[i][j][_n0-1][_n1][0] = TRparam[0].g_param[0]; // background first
							gfit_params_e_mpi[0].data[i][j][_n0-1][_n1][0] = TRparam[0].g_param[3*_n0+1]; // background error firstk

							for(_n2=1; _n2<4; _n2++)
							{
								gfit_params_mpi[0].data[i][j][_n0-1][_n1][_n2] = TRparam[0].g_param[3*_n1+_n2];
								//printf("%d %f\n", _n2, gfit_params_mpi[0].data[i][j][_n0-1][_n1][_n2]);
							}
							for(_n2=1; _n2<4; _n2++)
							{
								gfit_params_e_mpi[0].data[i][j][_n0-1][_n1][_n2] = TRparam[0].g_param[3*(_n0+_n1)+_n2+1];
							//	printf("%d %f\n", _n2, gfit_params_e_mpi[0].data[i][j][_n0-1][_n1][_n2]);
							}
						}

						// 1.1. calculate and save BIC
						gfit_BIC = -2.0*TRparam[0].maxLogLikeF + (3*_n0+1)*log(N_vel_data);
						if(gfit_BIC > 1E5) gfit_BIC = 1E3; // C integer overflow...
						//gfit_BIC = -2.0*TRparam[0].maxLogLikeF + (3*_n0+1)*(log(N_vel_data)-log(2*M_PI));
						//gfit_BIC = -2.0*TRparam[0].maxLogLikeF + 2*(3*_n0+1);
						gfit_params_mpi[0].data[i][j][_n0-1][0][4] = gfit_BIC; // save it to the array for the first Gaussian fit, distinguished by the number of Gaussians 
						//printf("ngauss:%d bic:%f\n", TRparam[0].n_gauss, gfit_params_mpi[0].data[i][j][_n0-1][0][4]);


						gfit_bic_temp[_n0-1][0] = (double)_n0; // # of gaussian functions fitted
						gfit_bic_temp[_n0-1][1] = gfit_BIC*1E3; // BIC : multiplied by 1E3 for using int comparison function

						// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
						// derive std of the single gaussian fit to decide whether to go multiple gaussian fit..
						// if s/n is less than 3 then do only single gaussian fit 
						if(_n0 == 1)
						{
							// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
							// derive std of residuals of the single gaussian fit: if s/n is less than 2 then the rest of the gaussian fits with the larger numnber of gaussians is not proceeded...
							vel_min = (0 - TRparam[0].velocity_offset)*TRparam[0].Velocity_Seperation + TRparam[0].Reference_Velocity;
							vel_max = (1*nax3-1 - TRparam[0].velocity_offset)*TRparam[0].Velocity_Seperation + TRparam[0].Reference_Velocity;
							if(vel_min > vel_max)
							{
								vel_temp = vel_min;
								vel_min = vel_max;
								vel_max = vel_temp;
							}

							// set the limits for computing the std based on a single Gaussian fit : this is for being used below
							// Xo is [0-1] in normalised values.
							Xo_lower_std = gfit_params_mpi[0].data[i][j][0][0][3] - 3*gfit_params_mpi[0].data[i][j][0][0][2];
							Xo_upper_std = gfit_params_mpi[0].data[i][j][0][0][3] + 3*gfit_params_mpi[0].data[i][j][0][0][2];

							// update std using the data in the outer regions..
							std = 0;
							std_n = 0;
							for(k=0; k<1*nax3; k++) // residuals for the original vel. section
							{
								vel = (k - TRparam[0].velocity_offset)*TRparam[0].Velocity_Seperation + TRparam[0].Reference_Velocity;
								vel_norm = (vel-vel_min)/(vel_max-vel_min); // normalise velocity to make it easier to fit

								//if(TRparam[0].Velocity_Seperation < 0) vel_norm = 1 - vel_norm;
								//gauss_model = gauss_function_at_a_pixel(TRparam[0].g_param, TRparam[0].n_gauss, _i_, _j_, vel_norm); // calculate a model value for a given polynomial function
								if((vel_norm < Xo_lower_std || vel_norm > Xo_upper_std))
								{
									std_n++;
									gauss_model = gauss_function_at_a_pixel_gfit_params_mpi(gfit_params_mpi[0].data, 1, i, j, vel_norm);

									velocity_profile_temp[std_n] = TRparam[0].velocity_profile_norm[k] - gauss_model;
								}
							}
							//std = sqrt(std / (std_n-1)); // standard deviation of residuals

							if(std_n < 3)
							{
								std = 0.1; // default value
							}
							else
							{
								robust_mean_std(velocity_profile_temp, std_n, &res_mean, &res_std);
								std = res_std;
							}

							// computing s/n
							_Xo = gfit_params_mpi[0].data[i][j][0][0][3];
							sn_profile = (each_gauss_function_at_a_pixel_gfit_params_mpi(gfit_params_mpi[0].data, 1, 0, i, j, _Xo) - gfit_params_mpi[0].data[i][j][0][0][0]) / std;
							gfit_params_mpi[0].data[i][j][0][0][7] = sn_profile;

							// put the std value
							gfit_params_mpi[0].data[i][j][0][0][6] = std;
							// +++++++++++++++++++++++++++++++++++++
							if(sn_profile < 0.5)
							{
								_n0 = n_gauss_upto + 2; // no more gaussian fit... exit loop
								n_gauss_upto_t = 1; // for sorting below
								gfit_params_mpi[0].data[i][j][0][0][4] = gfit_bic_temp[0][0];
								gfit_params_mpi[0].data[i][j][0][0][5] = 1; // single gaussian 
							}
						}
					}


					// save the optimal gauss-# based on BIC
					if(_n0 != (n_gauss_upto + 2))
					{
						qsort(gfit_bic_temp, n_gauss_upto_t, sizeof(gfit_bic_temp[0]), array2D_comp);
						gfit_params_mpi[0].data[i][j][0][0][5] = gfit_bic_temp[0][0];
						//gfit_params_mpi[0].data[i][j][0][0][5] = 4;
					}

					// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
					// derive std of residuals;
					vel_min = (0 - TRparam[0].velocity_offset)*TRparam[0].Velocity_Seperation + TRparam[0].Reference_Velocity;
					vel_max = (1*nax3-1 - TRparam[0].velocity_offset)*TRparam[0].Velocity_Seperation + TRparam[0].Reference_Velocity;
					if(vel_min > vel_max)
					{
						vel_temp = vel_min;
						vel_min = vel_max;
						vel_max = vel_temp;
					}

					std = 0;
					std_n = 0;
					for(k=0; k<1*nax3; k++) // residuals for the original vel. section
					{
						vel = (k - TRparam[0].velocity_offset)*TRparam[0].Velocity_Seperation + TRparam[0].Reference_Velocity;
						vel_norm = (vel-vel_min)/(vel_max-vel_min); // normalise velocity to make it easier to fit
						if((vel_norm < Xo_lower_std || vel_norm > Xo_upper_std))
						{
							std_n++;
							gauss_model = gauss_function_at_a_pixel_gfit_params_mpi(gfit_params_mpi[0].data, (int)gfit_params_mpi[0].data[i][j][0][0][5], i, j, vel_norm);
							velocity_profile_temp[std_n] = TRparam[0].velocity_profile_norm[k] - gauss_model;
						}

					}
					robust_mean_std(velocity_profile_temp, std_n, &res_mean, &res_std);
					std = res_std;

					if(std < 1E-4) // too small rms : this is for model cubes..  in unit 
					{
						std = 0;
						std_n = 0;
						for(k=0; k<1*nax3; k++) // residuals for the original vel. section
						{
							vel = (k - TRparam[0].velocity_offset)*TRparam[0].Velocity_Seperation + TRparam[0].Reference_Velocity;
							vel_norm = (vel-vel_min)/(vel_max-vel_min); // normalise velocity to make it easier to fit
							//if(TRparam[0].Velocity_Seperation < 0) vel_norm = 1 - vel_norm;
							//gauss_model = gauss_function_at_a_pixel(TRparam[0].g_param, TRparam[0].n_gauss, _i_, _j_, vel_norm); // calculate a model value for a given polynomial function
							if((k < (0+0.3*nax3) || k > (1*nax3-0.3*nax3))) // +- 10% of the edge regions..
							{
								std_n++;
								gauss_model = gauss_function_at_a_pixel_gfit_params_mpi(gfit_params_mpi[0].data, (int)gfit_params_mpi[0].data[i][j][0][0][5], i, j, vel_norm);

								velocity_profile_temp[std_n] = TRparam[0].velocity_profile_norm[k] - gauss_model;
							}
						}
						//std = sqrt(std / (std_n-1)); // standard deviation of residuals

						robust_mean_std(velocity_profile_temp, std_n, &res_mean, &res_std);
						std = res_std;
					}

					//for(k=0; k<nax3; k++) // update profile error with the std of the residuals between the observed and best derived gaussian fit results from the first fit
					//{
					//	TRparam[0].velocity_profile_e[i] = std;
						//printf("%d %d %f \n",i, j, std);
					//}
					// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

					// computing s/n for the profile whose amp is the lowest/highest among the deomposed the gaussian components..
					// +++++++++++++++++++++++++++++++++++++
					//_Xo = gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][0][3];
					for(_n3=0; _n3<(int)gfit_params_mpi[0].data[i][j][0][0][5]; _n3++)
					{
						gfit_A_temp[_n3][0] = gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][3]; // Xo of _n3-gauss 
						//gfit_A_temp[_n3][1] = 1E5*gauss_function_at_a_pixel_gfit_params_mpi(gfit_params_mpi[0].data, (int)gfit_params_mpi[0].data[i][j][0][0][5], i, j, gfit_A_temp[_n3][0]); // A at Xo : multiplied by 1E5 for using compare (int comparison function as it can't compare double values)
						gfit_A_temp[_n3][1] = 1E5*each_gauss_function_at_a_pixel_gfit_params_mpi(gfit_params_mpi[0].data, (int)gfit_params_mpi[0].data[i][j][0][0][5], _n3, i, j, gfit_A_temp[_n3][0]); // A at Xo : multiplied by 1E5 for using compare (int comparison function as it can't compare double values)
					}
					// sort gfit_A_temp to find the biggest A value which is used for computing s/n
					qsort(gfit_A_temp, (int)gfit_params_mpi[0].data[i][j][0][0][5], sizeof(gfit_A_temp[0]), array2D_comp);
					// the lowest flux [9]
					sn_profile = (gfit_A_temp[0][1]/1E5 - gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][0][0]) / std;
					gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][0][9] = sn_profile; // S/N of the profile with the LOWEST flux

					// the highest flux [8] : into the first G???1 component!
					sn_profile = (gfit_A_temp[(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][1]/1E5 - gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][0][0]) / std;
					gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][0][8] = sn_profile; // S/N of the profile with the HIGHEST flux

					// put optimal # of gaussians [5] : into the first G???1 component!
					gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][0][5] = gfit_params_mpi[0].data[i][j][0][0][5];

					// put the std value [6] : converted to the original units : into the first G???1 component!
					gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][0][6] = std*(amp_profile_max-amp_profile_min);

					// s/n of each gaussian : see below : [7]
					// +++++++++++++++++++++++++++++++++++++

					// computing s/n for each gaussian components : into all the G???? components
					// +++++++++++++++++++++++++++++++++++++
					for(_n3=0; _n3<(int)gfit_params_mpi[0].data[i][j][0][0][5]; _n3++)
					{
						_Xo = gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][3];
						sn_profile = (each_gauss_function_at_a_pixel_gfit_params_mpi(gfit_params_mpi[0].data, (int)gfit_params_mpi[0].data[i][j][0][0][5], _n3, i, j, _Xo) - gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][0][0]) / std;
						gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][7] = sn_profile;
					}
					// +++++++++++++++++++++++++++++++++++++
				}

				TRparam[0].n_gauss = n_gauss_upto; // recovered to the original one 
				n_gauss_upto_t = n_gauss_upto; // recovered to the input one

				// re-scaled peak flux at Xo : a which is converted to the Gaussian A below by multiplying sqrt(2*pi)*sigma where sigma is also converted to the original units: see below
				// here a is the peak flux of the each gaussian components in the original units conveted..	
				// below a is converted to the A : area under the each Gaussian comps. which can be converted to the integrated intensity map: moment 0
				//for(_n3=0; _n3<(int)gfit_params_mpi[0].data[i][j][0][0][5]; _n3++)
				//{
				//	_Xo = gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][3];
				//	gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][1] = each_gauss_function_at_a_pixel_gfit_params_mpi(gfit_params_mpi[0].data, (int)gfit_params_mpi[0].data[i][j][0][0][5], _n3, i, j, _Xo) * (amp_profile_max - amp_profile_min) + amp_profile_min;
				//	gfit_params_e_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][1] = gfit_params_e_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][1] * (amp_profile_max - amp_profile_min);
//printf("%d %d %//f %f %f\n", i, j, gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][1], gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][2], gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][3]);
				//}

				// +++++++++++++++++++++++++
				// recover the normalised units into the original ones... !!! this is for optimally decomposed gaussians found from Gfit_multinest!!!
				for(_n3=0; _n3<(int)gfit_params_mpi[0].data[i][j][0][0][5]; _n3++)
				{
					// background : to original units
					g0_norm = gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][0];
					gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][0] = gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][0] * (amp_profile_max - amp_profile_min) + amp_profile_min; // background
					gfit_params_e_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][0] = gfit_params_e_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][0] * (amp_profile_max - amp_profile_min); // background error

					// sigma
					s_norm = gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][2]; // save for the conversion of A below
					gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][2] = gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][2] * (vel_max - vel_min);
					gfit_params_e_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][2] = gfit_params_e_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][2] * (vel_max - vel_min);

					// intensity (gaussian area = A parameter)
					a_norm = gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][1];
					gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][1] = ((a_norm/(sqrt(2*M_PI)*s_norm) + g0_norm) * (amp_profile_max - amp_profile_min) + amp_profile_min -gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][0]) * sqrt(2*M_PI)*gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][2];
					//printf(" %d %d %f %f %f %f %f %d\n", i, j, amp_profile_max, amp_profile_min, gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][1], a_norm, (int)gfit_params_mpi[0].data[i][j][0][0][5]);
					ae_norm = gfit_params_e_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][1];
					gfit_params_e_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][1] = ((ae_norm/(sqrt(2*M_PI)*s_norm) + g0_norm) * (amp_profile_max - amp_profile_min) + amp_profile_min -gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][0]) * sqrt(2*M_PI)*gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][2];

					//gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][1] = gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][1] * sqrt(2.0*M_PI) * gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][2];
					//gfit_params_e_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][1] = gfit_params_e_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][1] * sqrt(2.0*M_PI) * gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][2];

					// central velocity
					gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][3] = gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][3] * (vel_max - vel_min) + vel_min;
					gfit_params_e_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][3] = gfit_params_e_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][3] * (vel_max - vel_min);


//printf("%d %d bg:%f a:%f s:%f v:%f %f\n", i, j, gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][0],
//									gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][1],
//									gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][2],
//									gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][_n3][3], a_norm);
				}
				// put the std value [6]
				//gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][0][6] = gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][0][6] * (amp_profile_max - amp_profile_min);
				// +++++++++++++++++++++++++

				// +++++++++++++++++++++++++
				// recover the normalised units into the original ones... !!! this is for the single gaussian component from Gfit_multinest!!! gfit_params_mpi[0].data[i][j][0][0][0-9]
				if((int)gfit_params_mpi[0].data[i][j][0][0][5] > 1) // if not done above.. : into the single component : sgfit...
				{
					// background : to original units
					g0_norm = gfit_params_mpi[0].data[i][j][0][0][0];
					gfit_params_mpi[0].data[i][j][0][0][0] = gfit_params_mpi[0].data[i][j][0][0][0] * (amp_profile_max - amp_profile_min) + amp_profile_min; // background
					gfit_params_e_mpi[0].data[i][j][0][0][0] = gfit_params_e_mpi[0].data[i][j][0][0][0] * (amp_profile_max - amp_profile_min); // background error

					// sigma
					s_norm = gfit_params_mpi[0].data[i][j][0][0][2];
					gfit_params_mpi[0].data[i][j][0][0][2] = gfit_params_mpi[0].data[i][j][0][0][2] * (vel_max - vel_min);
					gfit_params_e_mpi[0].data[i][j][0][0][2] = gfit_params_e_mpi[0].data[i][j][0][0][2] * (vel_max - vel_min);

					// intensity (gaussian area = A parameter) : A = a * sqrt(2pi) * sigma where sigma is in the original units
				
					a_norm = gfit_params_mpi[0].data[i][j][0][0][1];

					gfit_params_mpi[0].data[i][j][0][0][1] = ((a_norm/(sqrt(2*M_PI)*s_norm) + g0_norm) * (amp_profile_max - amp_profile_min) + amp_profile_min -gfit_params_mpi[0].data[i][j][0][0][0]) * sqrt(2*M_PI)*gfit_params_mpi[0].data[i][j][0][0][2]; // jy * m/s


					ae_norm = gfit_params_e_mpi[0].data[i][j][0][0][1];
					//gfit_params_mpi[0].data[i][j][0][0][1] = gfit_params_mpi[0].data[i][j][0][0][1] * sqrt(2.0*M_PI) * gfit_params_mpi[0].data[i][j][0][0][2];
					gfit_params_e_mpi[0].data[i][j][0][0][1] = ((ae_norm/(sqrt(2*M_PI)*s_norm) + g0_norm) * (amp_profile_max - amp_profile_min) + amp_profile_min -gfit_params_mpi[0].data[i][j][0][0][0]) * sqrt(2*M_PI)*gfit_params_mpi[0].data[i][j][0][0][2];

					// central velocity
					gfit_params_mpi[0].data[i][j][0][0][3] = gfit_params_mpi[0].data[i][j][0][0][3] * (vel_max - vel_min) + vel_min;
					gfit_params_e_mpi[0].data[i][j][0][0][3] = gfit_params_e_mpi[0].data[i][j][0][0][3] * (vel_max - vel_min);

					// rms
					gfit_params_mpi[0].data[i][j][0][0][6] =  gfit_params_mpi[0].data[i][j][0][0][6] * (amp_profile_max - amp_profile_min);
				}
				// +++++++++++++++++++++++++

				//printf("\n %d %d LOW:%f HIGH:%f %d std:%f %f %f\n", i, j, gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][0][9], gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][0][8], (int)gfit_params_mpi[0].data[i][j][0][0][5], gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][0][6], \
gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][0][1]/(sqrt(2*M_PI)*gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][0][2]), \
gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][1][1]/(sqrt(2*M_PI)*gfit_params_mpi[0].data[i][j][(int)gfit_params_mpi[0].data[i][j][0][0][5]-1][1][2]));
			}
			else
			{
				// recover the normalised units into the original ones...
				for(_n0=1; _n0<n_gauss_upto+1; _n0++)
				{

				//printf("after: amp_max:%f amp_min:%f\n", amp_profile_max, amp_profile_min);
					for(_n1=0; _n1<_n0; _n1++)
					{
						gfit_params_mpi[0].data[i][j][_n0-1][_n1][0] = 1E90;
						gfit_params_e_mpi[0].data[i][j][_n0-1][_n1][0] = 1E90;

						for(_n2=1; _n2<4; _n2++)
						{
							if(_n2 == 1) // amplitude
							{
								//printf("xxxx %f %f %f %f\n", gfit_params_mpi[0].data[i][j][_n0-1][_n1][_n2], amp_profile_max, amp_profile_min, gfit_params_mpi[0].data[i][j][_n0-1][_n1][_n2] * (amp_profile_max - amp_profile_min) + amp_profile_min);
								gfit_params_mpi[0].data[i][j][_n0-1][_n1][_n2] = 1E90;
							}
							else if(_n2 == 2)// sigma
							{
								gfit_params_mpi[0].data[i][j][_n0-1][_n1][_n2] = 1E90;
								gfit_params_mpi[0].data[i][j][_n0-1][_n1][1] = 1E90;
							}
							else if(_n2 == 3)// velocity
							{
								gfit_params_mpi[0].data[i][j][_n0-1][_n1][_n2] = 1E90;
							}
						//	printf("%d %f\n", _n2, gfit_params_mpi[0].data[i][j][_n0-1][_n1][_n2]);
						}
						for(_n2=1; _n2<4; _n2++)
						{
							if(_n2 == 1) // amplitude error
							{
								gfit_params_e_mpi[0].data[i][j][_n0-1][_n1][_n2] = 1E90;
							}
							else if(_n2 == 2)// sigma
							{
								gfit_params_e_mpi[0].data[i][j][_n0-1][_n1][_n2] = 1E90;
								gfit_params_e_mpi[0].data[i][j][_n0-1][_n1][1] = 1E90;
							}
							else if(_n2 == 3)// velocity
							{
								gfit_params_e_mpi[0].data[i][j][_n0-1][_n1][_n2] = 1E90;
							}
						//	printf("%d %f\n", _n2, gfit_params_e_mpi[0].data[i][j][_n0-1][_n1][_n2]);
						}
					}
				}
			}
			//printf("\r rank-%2d : x[%4d] y[%4d] : %4.2f : (%10d/%10d)", rank, i, j, 100*n_processed_profiles/n_total_profiles, (int)n_processed_profiles, (int)n_total_profiles);
			//fflush(stdout);

			//j += TRparam[0].decimY-1; 
			//if(j == nax2 || j > nax2)
			//	j = nax2-2;
		}
		//i += TRparam[0].decimX-1; 
		//if(i == nax1 || i > nax1)
		//	i = nax1-2;
	}

	free(gfit_bic_temp[0]);
	free(gfit_bic_temp);
	free(gfit_A_temp[0]);
	free(gfit_A_temp);
}

int cmp(const void *vp, const void *vq)
{
    const double *p = vp;
    const double *q = vq;
    double diff = *p - *q;
    return ((diff >= 0.0) ? ((diff > 0.0) ? -1 : 0) : +1);
}




// --- End of line
