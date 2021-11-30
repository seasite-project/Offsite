#pragma once

#define VARIANT_ID 123456

#include <math.h>
#include "RHS_InverterChain.h"
#include "ODE_radauIIA7.h"
#include "DS_J.h"
#ifdef INSTRUMENT
#include "timesnap.h"
#endif

void
timestep (const int me, const int first, const int last, double t, double h)
{
#ifdef INSTRUMENT
  if (me == 0)
    {
      printf ("\n#ImplVariant-123456\n");
    }
#endif
#pragma omp barrier
//RHS_predictor %11
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
	 imin (B_RHS_predictor_lj, last - first + 1), limit_RHS_predictor_lj =
	 imax (first, last + 1 - B_RHS_predictor_lj); jj <= last;
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
	printf ("#Kernel=11\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=11\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
#pragma omp barrier
  for (int k = 0; k < 6; ++k)
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
	     imin (B_LC_RHS_lji, last - first + 1), limit_LC_RHS_lji =
	     imax (first, last + 1 - B_LC_RHS_lji); jj <= last;
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
//ApproxUpdate %5
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
	 imin (B_ApproxUpdate_ij, last - first + 1), limit_ApproxUpdate_ij =
	 imax (first, last + 1 - B_ApproxUpdate_ij); jj <= last;
	 bs_ApproxUpdate_ij = last + 1 - jj, limit_ApproxUpdate_ij = last)
      {
	for (; jj <= limit_ApproxUpdate_ij; jj += bs_ApproxUpdate_ij)
	  {
#pragma nounroll_and_jam
	    for (int i = 0; i < 4; ++i)
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
	printf ("#Kernel=5\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=5\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
}
