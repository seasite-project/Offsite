#pragma once

#define VARIANT_ID 263

#include <math.h>
#include "RHS_InverterChain.h"
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
      printf ("\n#ImplVariant-263\n");
    }
#endif
#pragma omp barrier
//RHS_predictor %10
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
	printf ("#Kernel=10\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=10\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
#pragma omp barrier
//LC %18
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
	tmp = 0.112999479323160 * F[0][j];
	tmp += -0.0403092207235200 * F[1][j];
	tmp += 0.0258023774203400 * F[2][j];
	tmp += -0.00990467650730000 * F[3][j];
	Y[0][j] = tmp * h + y[j];
	tmp = 0.234383995747400 * F[0][j];
	tmp += 0.206892573935360 * F[1][j];
	tmp += -0.0478571280485400 * F[2][j];
	tmp += 0.0160474228065200 * F[3][j];
	Y[1][j] = tmp * h + y[j];
	tmp = 0.216681784623250 * F[0][j];
	tmp += 0.406123263867370 * F[1][j];
	tmp += 0.189036518170060 * F[2][j];
	tmp += -0.0241821048998300 * F[3][j];
	Y[2][j] = tmp * h + y[j];
	tmp = 0.220462211176770 * F[0][j];
	tmp += 0.388193468843170 * F[1][j];
	tmp += 0.328844319980060 * F[2][j];
	tmp += 0.0625000000000000 * F[3][j];
	Y[3][j] = tmp * h + y[j];
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
#pragma omp barrier
//LC %18
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
	    tmp = 0.112999479323160 * F[0][j];
	    tmp += -0.0403092207235200 * F[1][j];
	    tmp += 0.0258023774203400 * F[2][j];
	    tmp += -0.00990467650730000 * F[3][j];
	    Y[0][j] = tmp * h + y[j];
	    tmp = 0.234383995747400 * F[0][j];
	    tmp += 0.206892573935360 * F[1][j];
	    tmp += -0.0478571280485400 * F[2][j];
	    tmp += 0.0160474228065200 * F[3][j];
	    Y[1][j] = tmp * h + y[j];
	    tmp = 0.216681784623250 * F[0][j];
	    tmp += 0.406123263867370 * F[1][j];
	    tmp += 0.189036518170060 * F[2][j];
	    tmp += -0.0241821048998300 * F[3][j];
	    Y[2][j] = tmp * h + y[j];
	    tmp = 0.220462211176770 * F[0][j];
	    tmp += 0.388193468843170 * F[1][j];
	    tmp += 0.328844319980060 * F[2][j];
	    tmp += 0.0625000000000000 * F[3][j];
	    Y[3][j] = tmp * h + y[j];
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
//RHS_Approx %5
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
      }
    for (int i = 1; i < 4; ++i)
      {
#pragma ivdep
	for (int j = first; j <= last; ++j)
	  {
	    dy[j] += b[i] * (eval_component (j, t + c[i] * h, Y[i]));
	  }
      }
#ifdef INSTRUMENT
#pragma omp barrier
    if (me == 0)
      {
	double T = time_snap_stop (&time);
#ifdef _OPENMP
	printf ("#Kernel=5\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=5\t#Threads=1\t%.20e\n", T / 1e9 / n);
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
