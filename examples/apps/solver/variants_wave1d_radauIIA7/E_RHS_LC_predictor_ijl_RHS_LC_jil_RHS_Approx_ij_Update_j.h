#pragma once

#define VARIANT_ID 321

#include <math.h>
#include "RHS_Wave1D.h"
#include "ODE_radauIIA7.h"
#include "DS_E.h"
#ifdef INSTRUMENT
#include "timesnap.h"
#endif

void
timestep (const int me, const int first, const int last, double t, double h)
{
#ifdef INSTRUMENT
  if (me == 0)
    {
      printf ("\n#ImplVariant-321\n");
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
	double **swp_tmp = Yprev;
	Yprev = Ycur;
	Ycur = swp_tmp;
      }

//RHS_LC %28
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
	    tmp[0] = 0.112999479323160 * f;
	    tmp[1] = 0.234383995747400 * f;
	    tmp[2] = 0.216681784623250 * f;
	    tmp[3] = 0.220462211176770 * f;
	    f = eval_component (j, t + 0.409466864440740 * h, Yprev[1]);
	    tmp[0] += -0.0403092207235200 * f;
	    tmp[1] += 0.206892573935360 * f;
	    tmp[2] += 0.406123263867370 * f;
	    tmp[3] += 0.388193468843170 * f;
	    f = eval_component (j, t + 0.787659461760850 * h, Yprev[2]);
	    tmp[0] += 0.0258023774203400 * f;
	    tmp[1] += -0.0478571280485400 * f;
	    tmp[2] += 0.189036518170060 * f;
	    tmp[3] += 0.328844319980060 * f;
	    f = eval_component (j, t + 1.00000000000000 * h, Yprev[3]);
	    tmp[0] += -0.00990467650730000 * f;
	    tmp[1] += 0.0160474228065200 * f;
	    tmp[2] += -0.0241821048998300 * f;
	    tmp[3] += 0.0625000000000000 * f;
	    Ycur[0][j] = tmp[0] * h + y[j];
	    Ycur[1][j] = tmp[1] * h + y[j];
	    Ycur[2][j] = tmp[2] * h + y[j];
	    Ycur[3][j] = tmp[3] * h + y[j];
	  }
#ifdef INSTRUMENT
#pragma omp barrier
	if (me == 0)
	  {
	    double T = time_snap_stop (&time);
#ifdef _OPENMP
	    printf ("#Kernel=28\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=28\t#Threads=1\t%.20e\n", T / 1e9 / n);
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

//RHS_Approx %11
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
	  (eval_component (j, t + 0.0885879595126800 * h, Yprev[0]));
      }
    for (int i = 1; i < 4; ++i)
      {
#pragma ivdep
	for (int j = first; j <= last; ++j)
	  {
	    dy[j] += b[i] * (eval_component (j, t + c[i] * h, Yprev[i]));
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
