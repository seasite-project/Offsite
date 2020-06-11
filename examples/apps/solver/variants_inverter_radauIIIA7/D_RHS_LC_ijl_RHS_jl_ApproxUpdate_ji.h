#pragma once

#define VARIANT_ID 36

#include <math.h>
#include "RHS_InverterChain.h"
#include "ODE_radauIIA7.h"
#include "Y2sn_pFs.h"
#ifdef INSTRUMENT
#include "timesnap.h"
#endif

void
timestep (const int me, const int first, const int last, double t, double h)
{
#ifdef INSTRUMENT
  if (me == 0)
    {
    }
#endif
#pragma omp barrier
//RHS_LC %18
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
	f[0] = eval_component (j, t + 0.0885879595126800 * h, y);
	Ycur[0][j] = 0.112999479323160 * f[0];
	Ycur[1][j] = 0.234383995747400 * f[0];
	Ycur[2][j] = 0.216681784623250 * f[0];
	Ycur[3][j] = 0.220462211176770 * f[0];
      }
#pragma nounroll_and_jam
    for (int i = 1; i < 4; ++i)
      {
#pragma ivdep
	for (int j = first; j <= last; ++j)
	  {
	    f[0] = eval_component (j, t + c[i] * h, y);
	    Ycur[0][j] += A[0][i] * f[0];
	    Ycur[1][j] += A[1][i] * f[0];
	    Ycur[2][j] += A[2][i] * f[0];
	    Ycur[3][j] += A[3][i] * f[0];
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
	printf ("#Kernel=18\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=18\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
#pragma omp barrier
  for (int k = 0; k < 5; ++k)
    {

      double **tmp = Yprev;
      Yprev = Ycur;
      Ycur = tmp;

//RHS_LC %18
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
	    f[0] = eval_component (j, t + 0.0885879595126800 * h, Yprev[0]);
	    Ycur[0][j] = 0.112999479323160 * f[0];
	    Ycur[1][j] = 0.234383995747400 * f[0];
	    Ycur[2][j] = 0.216681784623250 * f[0];
	    Ycur[3][j] = 0.220462211176770 * f[0];
	  }
#pragma nounroll_and_jam
	for (int i = 1; i < 4; ++i)
	  {
#pragma ivdep
	    for (int j = first; j <= last; ++j)
	      {
		f[0] = eval_component (j, t + c[i] * h, Yprev[i]);
		Ycur[0][j] += A[0][i] * f[0];
		Ycur[1][j] += A[1][i] * f[0];
		Ycur[2][j] += A[2][i] * f[0];
		Ycur[3][j] += A[3][i] * f[0];
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
	    printf ("#Kernel=18\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=18\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
#pragma omp barrier
    }

  double **tmp = Yprev;
  Yprev = Ycur;
  Ycur = tmp;

//RHS %21
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
	F[0][j] = eval_component (j, t + 0.0885879595126800 * h, Yprev[0]);
	F[1][j] = eval_component (j, t + 0.409466864440740 * h, Yprev[1]);
	F[2][j] = eval_component (j, t + 0.787659461760850 * h, Yprev[2]);
	F[3][j] = eval_component (j, t + 1.00000000000000 * h, Yprev[3]);
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
//ApproxUpdate %4
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
	y[j] += h * 0.220462211176770 * F[0][j];
	y[j] += h * 0.388193468843170 * F[1][j];
	y[j] += h * 0.328844319980060 * F[2][j];
	y[j] += h * 0.0625000000000000 * F[3][j];
      }
#ifdef INSTRUMENT
#pragma omp barrier
    if (me == 0)
      {
	double T = time_snap_stop (&time);
#ifdef _OPENMP
	printf ("#Kernel=4\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=4\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
}