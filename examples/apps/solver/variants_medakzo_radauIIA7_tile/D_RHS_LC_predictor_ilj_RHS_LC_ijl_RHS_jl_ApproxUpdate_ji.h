#pragma once

#define VARIANT_ID 108

#include <math.h>
#include "RHS_Medakzo.h"
#include "ODE_radauIIA7.h"
#include "DS_D.h"
#ifdef INSTRUMENT
#include "timesnap.h"
#endif

void
timestep (const int me, const int first, const int last, double t, double h)
{
#ifdef INSTRUMENT
  if (me == 0)
    {
      printf ("\n#ImplVariant-108\n");
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

//RHS %33
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_RHS_jl, limit_RHS_jl;
    for (int jj = first, bs_RHS_jl =
	 imin (B, last - first + 1), limit_RHS_jl =
	 imax (first, last + 1 - B); jj <= last;
	 bs_RHS_jl = last + 1 - jj, limit_RHS_jl = last)
      {
	for (; jj <= limit_RHS_jl; jj += bs_RHS_jl)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_RHS_jl + jj; ++j)
	      {
		F[0][j] =
		  eval_component (j, t + 0.0885879595126800 * h, Yprev[0]);
		F[1][j] =
		  eval_component (j, t + 0.409466864440740 * h, Yprev[1]);
		F[2][j] =
		  eval_component (j, t + 0.787659461760850 * h, Yprev[2]);
		F[3][j] =
		  eval_component (j, t + 1.00000000000000 * h, Yprev[3]);
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
#pragma omp barrier
//ApproxUpdate %10
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_ApproxUpdate_ji, limit_ApproxUpdate_ji;
    for (int jj = first, bs_ApproxUpdate_ji =
	 imin (B, last - first + 1), limit_ApproxUpdate_ji =
	 imax (first, last + 1 - B); jj <= last;
	 bs_ApproxUpdate_ji = last + 1 - jj, limit_ApproxUpdate_ji = last)
      {
	for (; jj <= limit_ApproxUpdate_ji; jj += bs_ApproxUpdate_ji)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_ApproxUpdate_ji + jj; ++j)
	      {
		y[j] = h * 0.220462211176770 * F[0][j];
		y[j] += h * 0.388193468843170 * F[1][j];
		y[j] += h * 0.328844319980060 * F[2][j];
		y[j] += h * 0.0625000000000000 * F[3][j];
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
}
