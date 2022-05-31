#pragma once

#define VARIANT_ID 387

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
      printf ("\n#ImplVariant-387\n");
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
    for (int l = 0; l < 4; ++l)
      {
	eval_range (first, last, t + c[l] * h, y, F[l]);
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
#pragma nounroll_and_jam
	for (int l = 0; l < 4; ++l)
	  {
#pragma ivdep
	    for (int j = first; j <= last; ++j)
	      {
		Y[j] = A[l][0] * Fprev[0][j];
		Y[j] += A[l][1] * Fprev[1][j];
		Y[j] += A[l][2] * Fprev[2][j];
		Y[j] += A[l][3] * Fprev[3][j];
		Y[j] = Y[j] * h + y[j];
	      }
#pragma omp barrier
#pragma ivdep
	    for (int j = first; j <= last; ++j)
	      {
		F[l][j] = eval_component (j, t + c[l] * h, Y);
	      }
#pragma omp barrier
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
#pragma ivdep
    for (int j = first; j <= last; ++j)
      {
	y[j] = h * 0.220462211176770 * F[0][j];
      }
#pragma nounroll_and_jam
    for (int i = 1; i < 4; ++i)
      {
#pragma ivdep
	for (int j = first; j <= last; ++j)
	  {
	    y[j] += h * b[i] * F[i][j];
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
