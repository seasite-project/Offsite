#pragma once

#define VARIANT_ID 208

#include <math.h>
#include "RHS_Wave1D.h"
#include "ODE_radauIIA7.h"
#include "DS_C.h"
#ifdef INSTRUMENT
#include "timesnap.h"
#endif

void
timestep (const int me, const int first, const int last, double t, double h)
{
#ifdef INSTRUMENT
  if (me == 0)
    {
      printf ("\n#ImplVariant-208\n");
    }
#endif
#pragma omp barrier
//RHS_LC_predictor %14
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
	printf ("#Kernel=14\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=14\t#Threads=1\t%.20e\n", T / 1e9 / n);
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

//RHS_LC %25
#ifdef INSTRUMENT
      {
#pragma omp barrier
	time_snap_t time;
	if (me == 0)
	  {
	    time_snap_start (&time);
	  }
#endif
	int bs_RHS_LC_ilj, limit_RHS_LC_ilj;
	for (int jj = first, bs_RHS_LC_ilj =
	     imin (B, last - first + 1), limit_RHS_LC_ilj =
	     imax (first, last + 1 - B); jj <= last;
	     bs_RHS_LC_ilj = last + 1 - jj, limit_RHS_LC_ilj = last)
	  {
	    for (; jj <= limit_RHS_LC_ilj; jj += bs_RHS_LC_ilj)
	      {
		eval_range (first, last, t + 0.0885879595126800 * h, Yprev[0],
			    Fn);
		for (int l = 0; l < 4; ++l)
		  {
#pragma ivdep
		    for (int j = jj; j < bs_RHS_LC_ilj + jj; ++j)
		      {
			Ycur[l][j] = A[l][0] * Fn[j];
		      }
		  }
#pragma nounroll_and_jam
		for (int i = 1; i < 4; ++i)
		  {
		    eval_range (first, last, t + c[i] * h, Yprev[i], Fn);
		    for (int l = 0; l < 4; ++l)
		      {
#pragma ivdep
			for (int j = jj; j < bs_RHS_LC_ilj + jj; ++j)
			  {
			    Ycur[l][j] += A[l][i] * Fn[j];
			  }
		      }
		  }
#pragma nounroll_and_jam
		for (int l = 0; l < 4; ++l)
		  {
#pragma ivdep
		    for (int j = jj; j < bs_RHS_LC_ilj + jj; ++j)
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
	    printf ("#Kernel=25\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=25\t#Threads=1\t%.20e\n", T / 1e9 / n);
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

//RHS %27
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
	printf ("#Kernel=27\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=27\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
#pragma omp barrier
//Approx %8
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_Approx_ji, limit_Approx_ji;
    for (int jj = first, bs_Approx_ji =
	 imin (B, last - first + 1), limit_Approx_ji =
	 imax (first, last + 1 - B); jj <= last;
	 bs_Approx_ji = last + 1 - jj, limit_Approx_ji = last)
      {
	for (; jj <= limit_Approx_ji; jj += bs_Approx_ji)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_Approx_ji + jj; ++j)
	      {
		dy[j] = 0.220462211176770 * F[0][j];
		dy[j] += 0.388193468843170 * F[1][j];
		dy[j] += 0.328844319980060 * F[2][j];
		dy[j] += 0.0625000000000000 * F[3][j];
	      }
	  }
      }
#ifdef INSTRUMENT
#pragma omp barrier
    if (me == 0)
      {
	double T = time_snap_stop (&time);
#ifdef _OPENMP
	printf ("#Kernel=8\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=8\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
//Update %21
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
	printf ("#Kernel=21\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=21\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
}
