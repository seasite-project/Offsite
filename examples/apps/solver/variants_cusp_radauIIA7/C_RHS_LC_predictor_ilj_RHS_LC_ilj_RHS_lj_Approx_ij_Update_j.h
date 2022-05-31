#pragma once

#define VARIANT_ID 205

#include <math.h>
#include "RHS_Cusp.h"
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
      printf ("\n#ImplVariant-205\n");
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
    eval_range (first, last, t + 0.0885879595126800 * h, y, Fn);
    for (int l = 0; l < 4; ++l)
      {
#pragma ivdep
	for (int j = first; j <= last; ++j)
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
	    for (int j = first; j <= last; ++j)
	      {
		Ycur[l][j] += A[l][i] * Fn[j];
	      }
	  }
      }
#pragma nounroll_and_jam
    for (int l = 0; l < 4; ++l)
      {
#pragma ivdep
	for (int j = first; j <= last; ++j)
	  {
	    Ycur[l][j] = Ycur[l][j] * h + y[j];
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
	eval_range (first, last, t + 0.0885879595126800 * h, Yprev[0], Fn);
	for (int l = 0; l < 4; ++l)
	  {
#pragma ivdep
	    for (int j = first; j <= last; ++j)
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
		for (int j = first; j <= last; ++j)
		  {
		    Ycur[l][j] += A[l][i] * Fn[j];
		  }
	      }
	  }
#pragma nounroll_and_jam
	for (int l = 0; l < 4; ++l)
	  {
#pragma ivdep
	    for (int j = first; j <= last; ++j)
	      {
		Ycur[l][j] = Ycur[l][j] * h + y[j];
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

//RHS %26
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
	eval_range (first, last, t + c[l] * h, Yprev[l], F[l]);
      }
#ifdef INSTRUMENT
#pragma omp barrier
    if (me == 0)
      {
	double T = time_snap_stop (&time);
#ifdef _OPENMP
	printf ("#Kernel=26\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=26\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
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
