#pragma once

#define VARIANT_ID 9

#include <math.h>
#include "RHS_Medakzo.h"
#include "ODE_radauIIA7.h"
#include "DS_K.h"
#ifdef INSTRUMENT
#include "timesnap.h"
#endif

void
timestep (const int me, const int first, const int last, double t, double h)
{
#ifdef INSTRUMENT
  if (me == 0)
    {
      printf ("\n#ImplVariant-9\n");
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

//LC_RHS_Approx %1
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
#pragma ivdep
	for (int j = first; j <= last; ++j)
	  {
	    dy[j] += b[l] * (eval_component (j, t + c[l] * h, Y));
	  }
      }
#ifdef INSTRUMENT
#pragma omp barrier
    if (me == 0)
      {
	double T = time_snap_stop (&time);
#ifdef _OPENMP
	printf ("#Kernel=1\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=1\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
#pragma omp barrier
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
#pragma nounroll_and_jam
#pragma ivdep
    for (int j = first; j <= last; ++j)
      {
	y[j] += h * dy[j];
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
