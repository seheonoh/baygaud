
#include "baygaud.gsl.h"

// 2DBAT user defined functions
// GSL related

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void robust_mean_std(double *input, int n, double *robust_mean, double *robust_std)
{
    int i, n_filtered=0;
    int n_loop=0;
    double robust_mean_t;
    double robust_std_t;
    double filter_l, filter_u, sum_sqrs;

    // 1. Derive individual errors : robust mode values based on histograms
    if(n == 0)
    {
        //printf("The number of input array is %d\n", n);
        //printf("Put 1E9 for now\n");
        robust_mean_t = 1E9;
        robust_std_t = 1E9;
        return;
    }
    else if(n == 1)
    {
        //printf("The number of input array is %d\n", n);
        robust_mean_t = input[0];
        robust_std_t = 1E9;
        return;
    }
    gsl_sort(input, 1, n);
    robust_mean_t = gsl_stats_median_from_sorted_data(input, 1, n);
    robust_std_t = gsl_stats_sd_m(input, 1, n, robust_mean_t);

    // first pass with 3 sigma level 
    filter_l = robust_mean_t - 3.0*robust_std_t;
    filter_u = robust_mean_t + 3.0*robust_std_t;
    if(filter_l >= filter_u)
    {
        filter_l = filter_l - 1E9;
        filter_u = filter_l + 1E9;  
    }

    for(i=0; i<n; i++)
    {
        if(input[i] < filter_l || input[i] > filter_u)
        {
            input[i] = robust_mean_t;
        }
    }

    gsl_sort(input, 1, n);
    robust_mean_t = gsl_stats_median_from_sorted_data(input, 1, n);
    robust_std_t = gsl_stats_sd_m(input, 1, n, robust_mean_t);
    n_loop = 0;
    while(1)
    {
        n_loop++;
        n_filtered=0;
        filter_l = robust_mean_t - 3*robust_std_t;
        filter_u = robust_mean_t + 3*robust_std_t;

        for(i=0; i<n; i++)
        {
            if(input[i] < filter_l || input[i] > filter_u)
            {
                n_filtered++;
                input[i] = robust_mean_t;
            }
        }
        if(n_filtered == 0)
            break;

        if(n_loop > 1E3)
            break;

        gsl_sort(input, 1, n);
        robust_mean_t = gsl_stats_median_from_sorted_data(input, 1, n);
        robust_std_t = gsl_stats_sd_m(input, 1, n, robust_mean_t);
    }

    if(robust_std_t == 0)
        robust_std_t = 1E9;

    *robust_mean = robust_mean_t;
    *robust_std = robust_std_t;
    //printf("mean_input: %f std_input: %f\n", robust_mean_t, robust_std_t);
    return;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// --- End of line




