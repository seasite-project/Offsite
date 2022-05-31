#pragma once

#define VARIANT_ID 73

#include <math.h>
#include "RHS_Wave1D.h"
#include "ODE_radauIIA7.h"
#include "DS_B.h"
#ifdef INSTRUMENT
#include "timesnap.h"
#endif

void
timestep (const int me, const int first, const int last, double t, double h)
{
#ifdef INSTRUMENT
  if (me == 0)
    {
      printf ("\n#ImplVariant-73\n");
    }
#endif
#pragma omp barrier
//RHS_predictor %10
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_RHS_predictor_jl, limit_RHS_predictor_jl;
    for (int jj = first, bs_RHS_predictor_jl =
	 imin (B, last - first + 1), limit_RHS_predictor_jl =
	 imax (first, last + 1 - B); jj <= last;
	 bs_RHS_predictor_jl = last + 1 - jj, limit_RHS_predictor_jl = last)
      {
	for (; jj <= limit_RHS_predictor_jl; jj += bs_RHS_predictor_jl)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_RHS_predictor_jl + jj; ++j)
	      {
		F[0][j] = eval_component (j, t + 0.0885879595126800 * h, y);
		F[1][j] = eval_component (j, t + 0.409466864440740 * h, y);
		F[2][j] = eval_component (j, t + 0.787659461760850 * h, y);
		F[3][j] = eval_component (j, t + 1.00000000000000 * h, y);
	      }
	  }
      }
#ifdef INSTRUMENT
#pragma omp barrier
    if (me == 0)
      {
	double T = time_snap_stop (&time);
#ifdef _OPENMP
	printf ("#Kernel=10\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=10\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
  for (int k = 0; k < 6; ++k)
    {
#pragma omp barrier
//LC %19
#ifdef INSTRUMENT
      {
#pragma omp barrier
	time_snap_t time;
	if (me == 0)
	  {
	    time_snap_start (&time);
	  }
#endif
	int bs_LC_lij, limit_LC_lij;
	for (int jj = first, bs_LC_lij =
	     imin (B, last - first + 1), limit_LC_lij =
	     imax (first, last + 1 - B); jj <= last;
	     bs_LC_lij = last + 1 - jj, limit_LC_lij = last)
	  {
	    for (; jj <= limit_LC_lij; jj += bs_LC_lij)
	      {
#pragma nounroll_and_jam
		for (int l = 0; l < 4; ++l)
		  {
#pragma ivdep
		    for (int j = jj; j < bs_LC_lij + jj; ++j)
		      {
			Y[l][j] = A[l][0] * F[0][j];
		      }
#pragma nounroll_and_jam
		    for (int i = 1; i < 4; ++i)
		      {
#pragma ivdep
			for (int j = jj; j < bs_LC_lij + jj; ++j)
			  {
			    Y[l][j] += A[l][i] * F[i][j];
			  }
		      }
#pragma ivdep
		    for (int j = jj; j < bs_LC_lij + jj; ++j)
		      {
			Y[l][j] = Y[l][j] * h + y[j];
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
	    printf ("#Kernel=19\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=19\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
#pragma omp barrier
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
    }
#pragma omp barrier
//ApproxUpdate %3
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_ApproxUpdate_ij, limit_ApproxUpdate_ij;
    for (int jj = first, bs_ApproxUpdate_ij =
	 imin (B, last - first + 1), limit_ApproxUpdate_ij =
	 imax (first, last + 1 - B); jj <= last;
	 bs_ApproxUpdate_ij = last + 1 - jj, limit_ApproxUpdate_ij = last)
      {
	for (; jj <= limit_ApproxUpdate_ij; jj += bs_ApproxUpdate_ij)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_ApproxUpdate_ij + jj; ++j)
	      {
		y[j] = h * 0.220462211176770 * F[0][j];
	      }
#pragma nounroll_and_jam
	    for (int i = 1; i < 4; ++i)
	      {
#pragma ivdep
		for (int j = jj; j < bs_ApproxUpdate_ij + jj; ++j)
		  {
		    y[j] += h * b[i] * F[i][j];
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
	printf ("#Kernel=3\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=3\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
}
