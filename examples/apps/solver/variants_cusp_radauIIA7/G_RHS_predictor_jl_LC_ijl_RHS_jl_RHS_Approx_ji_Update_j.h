#pragma once

#define VARIANT_ID 260

#include <math.h>
#include "RHS_Cusp.h"
#include "ODE_radauIIA7.h"
#include "DS_G.h"
#ifdef INSTRUMENT
#include "timesnap.h"
#endif

void
timestep (const int me, const int first, const int last, double t, double h)
{
#ifdef INSTRUMENT
  if (me == 0)
    {
      printf ("\n#ImplVariant-260\n");
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
//LC %21
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
	Y[0][j] = 0.112999479323160 * F[0][j];
	Y[1][j] = 0.234383995747400 * F[0][j];
	Y[2][j] = 0.216681784623250 * F[0][j];
	Y[3][j] = 0.220462211176770 * F[0][j];
      }
#pragma nounroll_and_jam
    for (int i = 1; i < 4; ++i)
      {
#pragma ivdep
	for (int j = first; j <= last; ++j)
	  {
	    Y[0][j] += A[0][i] * F[i][j];
	    Y[1][j] += A[1][i] * F[i][j];
	    Y[2][j] += A[2][i] * F[i][j];
	    Y[3][j] += A[3][i] * F[i][j];
	  }
      }
#pragma ivdep
    for (int j = first; j <= last; ++j)
      {
	Y[0][j] = Y[0][j] * h + y[j];
	Y[1][j] = Y[1][j] * h + y[j];
	Y[2][j] = Y[2][j] * h + y[j];
	Y[3][j] = Y[3][j] * h + y[j];
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
#pragma omp barrier
  for (int k = 0; k < 5; ++k)
    {
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
	    printf ("#Kernel=33\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=33\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
#pragma omp barrier
//LC %21
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
	    Y[0][j] = 0.112999479323160 * F[0][j];
	    Y[1][j] = 0.234383995747400 * F[0][j];
	    Y[2][j] = 0.216681784623250 * F[0][j];
	    Y[3][j] = 0.220462211176770 * F[0][j];
	  }
#pragma nounroll_and_jam
	for (int i = 1; i < 4; ++i)
	  {
#pragma ivdep
	    for (int j = first; j <= last; ++j)
	      {
		Y[0][j] += A[0][i] * F[i][j];
		Y[1][j] += A[1][i] * F[i][j];
		Y[2][j] += A[2][i] * F[i][j];
		Y[3][j] += A[3][i] * F[i][j];
	      }
	  }
#pragma ivdep
	for (int j = first; j <= last; ++j)
	  {
	    Y[0][j] = Y[0][j] * h + y[j];
	    Y[1][j] = Y[1][j] * h + y[j];
	    Y[2][j] = Y[2][j] * h + y[j];
	    Y[3][j] = Y[3][j] * h + y[j];
	  }
#ifdef INSTRUMENT
#pragma omp barrier
	if (me == 0)
	  {
	    double T = time_snap_stop (&time);
#ifdef _OPENMP
	    printf ("#Kernel=21\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=21\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
#pragma omp barrier
    }
//RHS_Approx %10
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
	dy[j] =
	  0.220462211176770 *
	  (eval_component (j, t + 0.0885879595126800 * h, Y[0]));
	dy[j] +=
	  0.388193468843170 *
	  (eval_component (j, t + 0.409466864440740 * h, Y[1]));
	dy[j] +=
	  0.328844319980060 *
	  (eval_component (j, t + 0.787659461760850 * h, Y[2]));
	dy[j] +=
	  0.0625000000000000 *
	  (eval_component (j, t + 1.00000000000000 * h, Y[3]));
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
