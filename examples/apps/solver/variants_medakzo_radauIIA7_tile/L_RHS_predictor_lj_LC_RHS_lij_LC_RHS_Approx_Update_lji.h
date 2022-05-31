#pragma once

#define VARIANT_ID 402

#include <math.h>
#include "RHS_Medakzo.h"
#include "ODE_radauIIA7.h"
#include "DS_L.h"
#ifdef INSTRUMENT
#include "timesnap.h"
#endif

void
timestep (const int me, const int first, const int last, double t, double h)
{
#ifdef INSTRUMENT
  if (me == 0)
    {
      printf ("\n#ImplVariant-402\n");
    }
#endif
#pragma omp barrier
//RHS_predictor %15
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_RHS_predictor_lj, limit_RHS_predictor_lj;
    for (int jj = first, bs_RHS_predictor_lj =
	 imin (B, last - first + 1), limit_RHS_predictor_lj =
	 imax (first, last + 1 - B); jj <= last;
	 bs_RHS_predictor_lj = last + 1 - jj, limit_RHS_predictor_lj = last)
      {
	for (; jj <= limit_RHS_predictor_lj; jj += bs_RHS_predictor_lj)
	  {
	    for (int l = 0; l < 4; ++l)
	      {
		eval_range (first, last, t + c[l] * h, y, F[l]);
	      }
	  }
      }
#ifdef INSTRUMENT
#pragma omp barrier
    if (me == 0)
      {
	double T = time_snap_stop (&time);
#ifdef _OPENMP
	printf ("#Kernel=15\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=15\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
#pragma omp barrier
  for (int k = 0; k < 5; ++k)
    {

#pragma omp master
      {
	double **swp_tmp = Fprev;
	Fprev = F;
	F = swp_tmp;
      }

//LC_RHS %3
#ifdef INSTRUMENT
      {
#pragma omp barrier
	time_snap_t time;
	if (me == 0)
	  {
	    time_snap_start (&time);
	  }
#endif
	int bs_LC_RHS_lij, limit_LC_RHS_lij;
	for (int jj = first, bs_LC_RHS_lij =
	     imin (B, last - first + 1), limit_LC_RHS_lij =
	     imax (first, last + 1 - B); jj <= last;
	     bs_LC_RHS_lij = last + 1 - jj, limit_LC_RHS_lij = last)
	  {
	    for (; jj <= limit_LC_RHS_lij; jj += bs_LC_RHS_lij)
	      {
#pragma nounroll_and_jam
		for (int l = 0; l < 4; ++l)
		  {
#pragma ivdep
		    for (int j = jj; j < bs_LC_RHS_lij + jj; ++j)
		      {
			Y[j] = A[l][0] * Fprev[0][j];
		      }
		    for (int i = 1; i < 4; ++i)
		      {
#pragma ivdep
			for (int j = jj; j < bs_LC_RHS_lij + jj; ++j)
			  {
			    Y[j] += A[l][i] * Fprev[i][j];
			  }
		      }
#pragma ivdep
		    for (int j = jj; j < bs_LC_RHS_lij + jj; ++j)
		      {
			Y[j] = Y[j] * h + y[j];
		      }
#pragma omp barrier
		    eval_range (first, last, t + c[l] * h, Y, F[l]);
#pragma omp barrier
		  }
	      }
	  }
#ifdef INSTRUMENT
#pragma omp barrier
	if (me == 0)
	  {
	    double T = time_snap_stop (&time);
#ifdef _OPENMP
	    printf ("#Kernel=3\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		    T / 1e9 / n);
#else
	    printf ("#Kernel=3\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
    }

#pragma omp master
  {
    double **swp_tmp = Fprev;
    Fprev = F;
    F = swp_tmp;
  }

//LC_RHS_Approx_Update %6
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_LC_RHS_Approx_Update_lji, limit_LC_RHS_Approx_Update_lji;
    for (int jj = first, bs_LC_RHS_Approx_Update_lji =
	 imin (B, last - first + 1), limit_LC_RHS_Approx_Update_lji =
	 imax (first, last + 1 - B); jj <= last;
	 bs_LC_RHS_Approx_Update_lji =
	 last + 1 - jj, limit_LC_RHS_Approx_Update_lji = last)
      {
	for (; jj <= limit_LC_RHS_Approx_Update_lji;
	     jj += bs_LC_RHS_Approx_Update_lji)
	  {
#pragma nounroll_and_jam
	    for (int l = 0; l < 4; ++l)
	      {
#pragma ivdep
		for (int j = jj; j < bs_LC_RHS_Approx_Update_lji + jj; ++j)
		  {
		    Y[j] = A[l][0] * Fprev[0][j];
		    Y[j] += A[l][1] * Fprev[1][j];
		    Y[j] += A[l][2] * Fprev[2][j];
		    Y[j] += A[l][3] * Fprev[3][j];
		    Y[j] = Y[j] * h + y[j];
		  }
#pragma ivdep
		for (int j = jj; j < bs_LC_RHS_Approx_Update_lji + jj; ++j)
		  {
		    y[j] += h * b[l] * (eval_component (j, t + c[l] * h, Y));
		  }
	      }
	  }
      }
#ifdef INSTRUMENT
#pragma omp barrier
    if (me == 0)
      {
	double T = time_snap_stop (&time);
#ifdef _OPENMP
	printf ("#Kernel=6\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=6\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
}
