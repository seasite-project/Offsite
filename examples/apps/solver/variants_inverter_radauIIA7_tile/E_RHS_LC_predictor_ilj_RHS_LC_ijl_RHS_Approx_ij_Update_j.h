#pragma once

#define VARIANT_ID 333

#include <math.h>
#include "RHS_InverterChain.h"
#include "ODE_radauIIA7.h"
#include "DS_E.h"
#ifdef INSTRUMENT
#include "timesnap.h"
#endif

void
timestep (const int me, const int first, const int last, double t, double h)
{
#ifdef INSTRUMENT
  if (me == 0)
    {
      printf ("\n#ImplVariant-333\n");
    }
#endif
#pragma omp barrier
//RHS_LC_predictor %20
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_RHS_LC_predictor_ilj, limit_RHS_LC_predictor_ilj;
    for (int jj = first, bs_RHS_LC_predictor_ilj =
	 imin (B, last - first + 1), limit_RHS_LC_predictor_ilj =
	 imax (first, last + 1 - B); jj <= last;
	 bs_RHS_LC_predictor_ilj = last + 1 - jj, limit_RHS_LC_predictor_ilj =
	 last)
      {
	for (; jj <= limit_RHS_LC_predictor_ilj;
	     jj += bs_RHS_LC_predictor_ilj)
	  {
	    eval_range (first, last, t + 0.0885879595126800 * h, y, Fn);
	    for (int l = 0; l < 4; ++l)
	      {
#pragma ivdep
		for (int j = jj; j < bs_RHS_LC_predictor_ilj + jj; ++j)
		  {
		    Ycur[l][j] = A[l][0] * Fn[j];
		  }
	      }
#pragma nounroll_and_jam
	    for (int i = 1; i < 4; ++i)
	      {
		eval_range (first, last, t + c[i] * h, y, Fn);
		for (int l = 0; l < 4; ++l)
		  {
#pragma ivdep
		    for (int j = jj; j < bs_RHS_LC_predictor_ilj + jj; ++j)
		      {
			Ycur[l][j] += A[l][i] * Fn[j];
		      }
		  }
	      }
#pragma nounroll_and_jam
	    for (int l = 0; l < 4; ++l)
	      {
#pragma ivdep
		for (int j = jj; j < bs_RHS_LC_predictor_ilj + jj; ++j)
		  {
		    Ycur[l][j] = Ycur[l][j] * h + y[j];
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

#pragma omp master
      {
	double **swp_tmp = Yprev;
	Yprev = Ycur;
	Ycur = swp_tmp;
      }

//RHS_LC %30
#ifdef INSTRUMENT
      {
#pragma omp barrier
	time_snap_t time;
	if (me == 0)
	  {
	    time_snap_start (&time);
	  }
#endif
	int bs_RHS_LC_ijl, limit_RHS_LC_ijl;
	for (int jj = first, bs_RHS_LC_ijl =
	     imin (B, last - first + 1), limit_RHS_LC_ijl =
	     imax (first, last + 1 - B); jj <= last;
	     bs_RHS_LC_ijl = last + 1 - jj, limit_RHS_LC_ijl = last)
	  {
	    for (; jj <= limit_RHS_LC_ijl; jj += bs_RHS_LC_ijl)
	      {
#pragma ivdep
		for (int j = jj; j < bs_RHS_LC_ijl + jj; ++j)
		  {
		    f =
		      eval_component (j, t + 0.0885879595126800 * h,
				      Yprev[0]);
		    Ycur[0][j] = 0.112999479323160 * f;
		    Ycur[1][j] = 0.234383995747400 * f;
		    Ycur[2][j] = 0.216681784623250 * f;
		    Ycur[3][j] = 0.220462211176770 * f;
		  }
#pragma nounroll_and_jam
		for (int i = 1; i < 4; ++i)
		  {
#pragma ivdep
		    for (int j = jj; j < bs_RHS_LC_ijl + jj; ++j)
		      {
			f = eval_component (j, t + c[i] * h, Yprev[i]);
			Ycur[0][j] += A[0][i] * f;
			Ycur[1][j] += A[1][i] * f;
			Ycur[2][j] += A[2][i] * f;
			Ycur[3][j] += A[3][i] * f;
		      }
		  }
#pragma ivdep
		for (int j = jj; j < bs_RHS_LC_ijl + jj; ++j)
		  {
		    Ycur[0][j] = Ycur[0][j] * h + y[j];
		    Ycur[1][j] = Ycur[1][j] * h + y[j];
		    Ycur[2][j] = Ycur[2][j] * h + y[j];
		    Ycur[3][j] = Ycur[3][j] * h + y[j];
		  }
	      }
	  }
#ifdef INSTRUMENT
#pragma omp barrier
	if (me == 0)
	  {
	    double T = time_snap_stop (&time);
#ifdef _OPENMP
	    printf ("#Kernel=30\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=30\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
#pragma omp barrier
    }

#pragma omp master
  {
    double **swp_tmp = Yprev;
    Yprev = Ycur;
    Ycur = swp_tmp;
  }

//RHS_Approx %11
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_RHS_Approx_ij, limit_RHS_Approx_ij;
    for (int jj = first, bs_RHS_Approx_ij =
	 imin (B, last - first + 1), limit_RHS_Approx_ij =
	 imax (first, last + 1 - B); jj <= last;
	 bs_RHS_Approx_ij = last + 1 - jj, limit_RHS_Approx_ij = last)
      {
	for (; jj <= limit_RHS_Approx_ij; jj += bs_RHS_Approx_ij)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_RHS_Approx_ij + jj; ++j)
	      {
		dy[j] =
		  0.220462211176770 *
		  (eval_component (j, t + 0.0885879595126800 * h, Yprev[0]));
	      }
	    for (int i = 1; i < 4; ++i)
	      {
#pragma ivdep
		for (int j = jj; j < bs_RHS_Approx_ij + jj; ++j)
		  {
		    dy[j] +=
		      b[i] * (eval_component (j, t + c[i] * h, Yprev[i]));
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
	printf ("#Kernel=11\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=11\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
//Update %27
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_Update_j, limit_Update_j;
    for (int jj = first, bs_Update_j =
	 imin (B, last - first + 1), limit_Update_j =
	 imax (first, last + 1 - B); jj <= last;
	 bs_Update_j = last + 1 - jj, limit_Update_j = last)
      {
	for (; jj <= limit_Update_j; jj += bs_Update_j)
	  {
#pragma nounroll_and_jam
#pragma ivdep
	    for (int j = jj; j < bs_Update_j + jj; ++j)
	      {
		y[j] += h * dy[j];
	      }
	  }
      }
#ifdef INSTRUMENT
#pragma omp barrier
    if (me == 0)
      {
	double T = time_snap_stop (&time);
#ifdef _OPENMP
	printf ("#Kernel=27\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=27\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
}
