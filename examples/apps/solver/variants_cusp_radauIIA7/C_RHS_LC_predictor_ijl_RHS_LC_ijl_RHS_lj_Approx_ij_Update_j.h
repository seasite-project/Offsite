#pragma once

#define VARIANT_ID 161

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
      printf ("\n#ImplVariant-161\n");
    }
#endif
#pragma omp barrier
//RHS_LC_predictor %19
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
	f = eval_component (j, t + 0.0885879595126800 * h, y);
	Ycur[0][j] = 0.112999479323160 * f;
	Ycur[1][j] = 0.234383995747400 * f;
	Ycur[2][j] = 0.216681784623250 * f;
	Ycur[3][j] = 0.220462211176770 * f;
      }
#pragma nounroll_and_jam
    for (int i = 1; i < 4; ++i)
      {
#pragma ivdep
	for (int j = first; j <= last; ++j)
	  {
	    f = eval_component (j, t + c[i] * h, y);
	    Ycur[0][j] += A[0][i] * f;
	    Ycur[1][j] += A[1][i] * f;
	    Ycur[2][j] += A[2][i] * f;
	    Ycur[3][j] += A[3][i] * f;
	  }
      }
#pragma ivdep
    for (int j = first; j <= last; ++j)
      {
	Ycur[0][j] = Ycur[0][j] * h + y[j];
	Ycur[1][j] = Ycur[1][j] * h + y[j];
	Ycur[2][j] = Ycur[2][j] * h + y[j];
	Ycur[3][j] = Ycur[3][j] * h + y[j];
      }
#ifdef INSTRUMENT
#pragma omp barrier
    if (me == 0)
      {
	double T = time_snap_stop (&time);
#ifdef _OPENMP
	printf ("#Kernel=19\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=19\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
#pragma omp barrier
  for (int k = 0; k < 5; ++k)
    {

#pragma omp master
      {
	double **tmp = Yprev;
	Yprev = Ycur;
	Ycur = tmp;
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
#pragma ivdep
	for (int j = first; j <= last; ++j)
	  {
	    f = eval_component (j, t + 0.0885879595126800 * h, Yprev[0]);
	    Ycur[0][j] = 0.112999479323160 * f;
	    Ycur[1][j] = 0.234383995747400 * f;
	    Ycur[2][j] = 0.216681784623250 * f;
	    Ycur[3][j] = 0.220462211176770 * f;
	  }
#pragma nounroll_and_jam
	for (int i = 1; i < 4; ++i)
	  {
#pragma ivdep
	    for (int j = first; j <= last; ++j)
	      {
		f = eval_component (j, t + c[i] * h, Yprev[i]);
		Ycur[0][j] += A[0][i] * f;
		Ycur[1][j] += A[1][i] * f;
		Ycur[2][j] += A[2][i] * f;
		Ycur[3][j] += A[3][i] * f;
	      }
	  }
#pragma ivdep
	for (int j = first; j <= last; ++j)
	  {
	    Ycur[0][j] = Ycur[0][j] * h + y[j];
	    Ycur[1][j] = Ycur[1][j] * h + y[j];
	    Ycur[2][j] = Ycur[2][j] * h + y[j];
	    Ycur[3][j] = Ycur[3][j] * h + y[j];
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
    double **tmp = Yprev;
    Yprev = Ycur;
    Ycur = tmp;
  }

//RHS %32
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
	printf ("#Kernel=32\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=32\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
#pragma omp barrier
//Approx %13
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
	printf ("#Kernel=13\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=13\t#Threads=1\t%.20e\n", T / 1e9 / n);
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
