#pragma once

#define VARIANT_ID 372

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
      printf ("\n#ImplVariant-372\n");
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
	double **tmp = Fprev;
	Fprev = F;
	F = tmp;
      }

//LC_RHS %8
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
	    printf ("#Kernel=8\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		    T / 1e9 / n);
#else
	    printf ("#Kernel=8\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
    }

#pragma omp master
  {
    double **tmp = Fprev;
    Fprev = F;
    F = tmp;
  }

//LC_RHS_Approx_Update %12
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
	printf ("#Kernel=12\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=12\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
}
