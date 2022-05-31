#pragma once

#define VARIANT_ID 24

#include <math.h>
#include "RHS_Heat2D.h"
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
      printf ("\n#ImplVariant-24\n");
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
#pragma omp barrier
  for (int k = 0; k < 5; ++k)
    {

#pragma omp master
      {
	double **swp_tmp = Fprev;
	Fprev = F;
	F = swp_tmp;
      }

//LC_RHS %31
#ifdef INSTRUMENT
      {
#pragma omp barrier
	time_snap_t time;
	if (me == 0)
	  {
	    time_snap_start (&time);
	  }
#endif
	int bs_LC_RHS_lji, limit_LC_RHS_lji;
	for (int jj = first, bs_LC_RHS_lji =
	     imin (B, last - first + 1), limit_LC_RHS_lji =
	     imax (first, last + 1 - B); jj <= last;
	     bs_LC_RHS_lji = last + 1 - jj, limit_LC_RHS_lji = last)
	  {
	    for (; jj <= limit_LC_RHS_lji; jj += bs_LC_RHS_lji)
	      {
#pragma nounroll_and_jam
		for (int l = 0; l < 4; ++l)
		  {
#pragma ivdep
		    for (int j = jj; j < bs_LC_RHS_lji + jj; ++j)
		      {
			Y[j] = A[l][0] * Fprev[0][j];
			Y[j] += A[l][1] * Fprev[1][j];
			Y[j] += A[l][2] * Fprev[2][j];
			Y[j] += A[l][3] * Fprev[3][j];
			Y[j] = Y[j] * h + y[j];
		      }
#pragma omp barrier
#pragma ivdep
		    for (int j = jj; j < bs_LC_RHS_lji + jj; ++j)
		      {
			F[l][j] = eval_component (j, t + c[l] * h, Y);
		      }
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
	    printf ("#Kernel=31\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=31\t#Threads=1\t%.20e\n", T / 1e9 / n);
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

//LC_RHS_Approx_Update %33
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
	printf ("#Kernel=33\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=33\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
}
