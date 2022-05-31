#pragma once

#define VARIANT_ID 230

#include <math.h>
#include "RHS_Heat2D.h"
#include "ODE_radauIIA7.h"
#include "DS_H.h"
#ifdef INSTRUMENT
#include "timesnap.h"
#endif

void
timestep (const int me, const int first, const int last, double t, double h)
{
#ifdef INSTRUMENT
  if (me == 0)
    {
      printf ("\n#ImplVariant-230\n");
    }
#endif
#pragma omp barrier
//RHS_predictor %9
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
	printf ("#Kernel=9\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=9\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
#pragma omp barrier
//LC %20
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_LC_lji, limit_LC_lji;
    for (int jj = first, bs_LC_lji =
	 imin (B, last - first + 1), limit_LC_lji =
	 imax (first, last + 1 - B); jj <= last;
	 bs_LC_lji = last + 1 - jj, limit_LC_lji = last)
      {
	for (; jj <= limit_LC_lji; jj += bs_LC_lji)
	  {
#pragma nounroll_and_jam
	    for (int l = 0; l < 4; ++l)
	      {
#pragma ivdep
		for (int j = jj; j < bs_LC_lji + jj; ++j)
		  {
		    tmp = A[l][0] * F[0][j];
		    tmp += A[l][1] * F[1][j];
		    tmp += A[l][2] * F[2][j];
		    tmp += A[l][3] * F[3][j];
		    Y[l][j] = tmp * h + y[j];
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
	printf ("#Kernel=20\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=20\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
#pragma omp barrier
  for (int k = 0; k < 5; ++k)
    {
//RHS %26
#ifdef INSTRUMENT
      {
#pragma omp barrier
	time_snap_t time;
	if (me == 0)
	  {
	    time_snap_start (&time);
	  }
#endif
	int bs_RHS_lj, limit_RHS_lj;
	for (int jj = first, bs_RHS_lj =
	     imin (B, last - first + 1), limit_RHS_lj =
	     imax (first, last + 1 - B); jj <= last;
	     bs_RHS_lj = last + 1 - jj, limit_RHS_lj = last)
	  {
	    for (; jj <= limit_RHS_lj; jj += bs_RHS_lj)
	      {
		for (int l = 0; l < 4; ++l)
		  {
		    eval_range (first, last, t + c[l] * h, Y[l], F[l]);
		  }
	      }
	  }
#ifdef INSTRUMENT
#pragma omp barrier
	if (me == 0)
	  {
	    double T = time_snap_stop (&time);
#ifdef _OPENMP
	    printf ("#Kernel=26\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=26\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
#pragma omp barrier
//LC %20
#ifdef INSTRUMENT
      {
#pragma omp barrier
	time_snap_t time;
	if (me == 0)
	  {
	    time_snap_start (&time);
	  }
#endif
	int bs_LC_lji, limit_LC_lji;
	for (int jj = first, bs_LC_lji =
	     imin (B, last - first + 1), limit_LC_lji =
	     imax (first, last + 1 - B); jj <= last;
	     bs_LC_lji = last + 1 - jj, limit_LC_lji = last)
	  {
	    for (; jj <= limit_LC_lji; jj += bs_LC_lji)
	      {
#pragma nounroll_and_jam
		for (int l = 0; l < 4; ++l)
		  {
#pragma ivdep
		    for (int j = jj; j < bs_LC_lji + jj; ++j)
		      {
			tmp = A[l][0] * F[0][j];
			tmp += A[l][1] * F[1][j];
			tmp += A[l][2] * F[2][j];
			tmp += A[l][3] * F[3][j];
			Y[l][j] = tmp * h + y[j];
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
	    printf ("#Kernel=20\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=20\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
#pragma omp barrier
    }
//RHS_Approx_Update %2
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_RHS_Approx_Update_ij, limit_RHS_Approx_Update_ij;
    for (int jj = first, bs_RHS_Approx_Update_ij =
	 imin (B, last - first + 1), limit_RHS_Approx_Update_ij =
	 imax (first, last + 1 - B); jj <= last;
	 bs_RHS_Approx_Update_ij = last + 1 - jj, limit_RHS_Approx_Update_ij =
	 last)
      {
	for (; jj <= limit_RHS_Approx_Update_ij;
	     jj += bs_RHS_Approx_Update_ij)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_RHS_Approx_Update_ij + jj; ++j)
	      {
		y[j] =
		  h * 0.220462211176770 *
		  (eval_component (j, t + 0.0885879595126800 * h, Y[0]));
	      }
#pragma nounroll_and_jam
	    for (int i = 1; i < 4; ++i)
	      {
#pragma ivdep
		for (int j = jj; j < bs_RHS_Approx_Update_ij + jj; ++j)
		  {
		    y[j] +=
		      h * b[i] * (eval_component (j, t + c[i] * h, Y[i]));
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
	printf ("#Kernel=2\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=2\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
}
