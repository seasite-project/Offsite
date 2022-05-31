#pragma once

#define VARIANT_ID 22

#include <math.h>
#include "RHS_Wave1D.h"
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
      printf ("\n#ImplVariant-22\n");
    }
#endif
#pragma omp barrier
//RHS_predictor %16
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
	F[0][j] = eval_component (j, t + 0.0885879595126800 * h, y);
	F[1][j] = eval_component (j, t + 0.409466864440740 * h, y);
	F[2][j] = eval_component (j, t + 0.787659461760850 * h, y);
	F[3][j] = eval_component (j, t + 1.00000000000000 * h, y);
      }
#ifdef INSTRUMENT
#pragma omp barrier
    if (me == 0)
      {
	double T = time_snap_stop (&time);
#ifdef _OPENMP
	printf ("#Kernel=16\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=16\t#Threads=1\t%.20e\n", T / 1e9 / n);
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
#pragma nounroll_and_jam
	for (int l = 0; l < 4; ++l)
	  {
#pragma ivdep
	    for (int j = first; j <= last; ++j)
	      {
		Y[j] = A[l][0] * Fprev[0][j];
	      }
	    for (int i = 1; i < 4; ++i)
	      {
#pragma ivdep
		for (int j = first; j <= last; ++j)
		  {
		    Y[j] += A[l][i] * Fprev[i][j];
		  }
	      }
#pragma ivdep
	    for (int j = first; j <= last; ++j)
	      {
		Y[j] = Y[j] * h + y[j];
	      }
#pragma omp barrier
	    eval_range (first, last, t + c[l] * h, Y, F[l]);
#pragma omp barrier
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
#pragma ivdep
	for (int j = first; j <= last; ++j)
	  {
	    y[j] += h * b[l] * (eval_component (j, t + c[l] * h, Y));
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
