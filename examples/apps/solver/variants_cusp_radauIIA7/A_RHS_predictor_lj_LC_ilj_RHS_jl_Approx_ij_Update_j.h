#pragma once

#define VARIANT_ID 343

#include <math.h>
#include "RHS_Cusp.h"
#include "ODE_radauIIA7.h"
#include "DS_A.h"
#ifdef INSTRUMENT
#include "timesnap.h"
#endif

void
timestep (const int me, const int first, const int last, double t, double h)
{
#ifdef INSTRUMENT
  if (me == 0)
    {
      printf ("\n#ImplVariant-343\n");
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
  for (int k = 0; k < 6; ++k)
    {
#pragma omp barrier
//LC %16
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
		Y[l][j] = A[l][0] * F[0][j];
	      }
	  }
#pragma nounroll_and_jam
	for (int i = 1; i < 4; ++i)
	  {
#pragma nounroll_and_jam
	    for (int l = 0; l < 4; ++l)
	      {
#pragma ivdep
		for (int j = first; j <= last; ++j)
		  {
		    Y[l][j] += A[l][i] * F[i][j];
		  }
	      }
	  }
	for (int l = 0; l < 4; ++l)
	  {
#pragma ivdep
	    for (int j = first; j <= last; ++j)
	      {
		Y[l][j] = Y[l][j] * h + y[j];
	      }
	  }
#ifdef INSTRUMENT
#pragma omp barrier
	if (me == 0)
	  {
	    double T = time_snap_stop (&time);
#ifdef _OPENMP
	    printf ("#Kernel=16\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=16\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
#pragma omp barrier
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
#pragma ivdep
	for (int j = first; j <= last; ++j)
	  {
	    F[0][j] = eval_component (j, t + 0.0885879595126800 * h, Y[0]);
	    F[1][j] = eval_component (j, t + 0.409466864440740 * h, Y[1]);
	    F[2][j] = eval_component (j, t + 0.787659461760850 * h, Y[2]);
	    F[3][j] = eval_component (j, t + 1.00000000000000 * h, Y[3]);
	  }
#ifdef INSTRUMENT
#pragma omp barrier
	if (me == 0)
	  {
	    double T = time_snap_stop (&time);
#ifdef _OPENMP
	    printf ("#Kernel=27\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=27\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
    }
#pragma omp barrier
//Approx %7
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
	dy[j] = 0.220462211176770 * F[0][j];
      }
#pragma nounroll_and_jam
    for (int i = 1; i < 4; ++i)
      {
#pragma ivdep
	for (int j = first; j <= last; ++j)
	  {
	    dy[j] += b[i] * F[i][j];
	  }
      }
#ifdef INSTRUMENT
#pragma omp barrier
    if (me == 0)
      {
	double T = time_snap_stop (&time);
#ifdef _OPENMP
	printf ("#Kernel=7\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=7\t#Threads=1\t%.20e\n", T / 1e9 / n);
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
	printf ("#Kernel=21\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=21\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
}
