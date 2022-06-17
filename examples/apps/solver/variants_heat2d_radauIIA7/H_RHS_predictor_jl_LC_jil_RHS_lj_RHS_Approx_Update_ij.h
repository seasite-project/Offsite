#pragma once

#define VARIANT_ID 210

#include <math.h>
#include "RHS_Heat2D.h"
#include "ODE_radauIIA7.h"
#include "DS_H.h"
#ifdef INSTRUMENT
#include "timesnap.h"
#endif

void
timestep (const int me, const int first, const int last, double t, double h)
{
#ifdef INSTRUMENT
  if (me == 0)
    {
      printf ("\n#ImplVariant-210\n");
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
//LC %23
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
	Yn[0] = 0.112999479323160 * F[0][j];
	Yn[1] = 0.234383995747400 * F[0][j];
	Yn[2] = 0.216681784623250 * F[0][j];
	Yn[3] = 0.220462211176770 * F[0][j];
	Yn[0] += -0.0403092207235200 * F[1][j];
	Yn[1] += 0.206892573935360 * F[1][j];
	Yn[2] += 0.406123263867370 * F[1][j];
	Yn[3] += 0.388193468843170 * F[1][j];
	Yn[0] += 0.0258023774203400 * F[2][j];
	Yn[1] += -0.0478571280485400 * F[2][j];
	Yn[2] += 0.189036518170060 * F[2][j];
	Yn[3] += 0.328844319980060 * F[2][j];
	Yn[0] += -0.00990467650730000 * F[3][j];
	Yn[1] += 0.0160474228065200 * F[3][j];
	Yn[2] += -0.0241821048998300 * F[3][j];
	Yn[3] += 0.0625000000000000 * F[3][j];
	Y[0][j] = Yn[0] * h + y[j];
	Y[1][j] = Yn[1] * h + y[j];
	Y[2][j] = Yn[2] * h + y[j];
	Y[3][j] = Yn[3] * h + y[j];
      }
#ifdef INSTRUMENT
#pragma omp barrier
    if (me == 0)
      {
	double T = time_snap_stop (&time);
#ifdef _OPENMP
	printf ("#Kernel=23\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=23\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
#pragma omp barrier
  for (int k = 0; k < 5; ++k)
    {
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
	    eval_range (first, last, t + c[l] * h, Y[l], F[l]);
	  }
#ifdef INSTRUMENT
#pragma omp barrier
	if (me == 0)
	  {
	    double T = time_snap_stop (&time);
#ifdef _OPENMP
	    printf ("#Kernel=32\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=32\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
#pragma omp barrier
//LC %23
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
	    Yn[0] = 0.112999479323160 * F[0][j];
	    Yn[1] = 0.234383995747400 * F[0][j];
	    Yn[2] = 0.216681784623250 * F[0][j];
	    Yn[3] = 0.220462211176770 * F[0][j];
	    Yn[0] += -0.0403092207235200 * F[1][j];
	    Yn[1] += 0.206892573935360 * F[1][j];
	    Yn[2] += 0.406123263867370 * F[1][j];
	    Yn[3] += 0.388193468843170 * F[1][j];
	    Yn[0] += 0.0258023774203400 * F[2][j];
	    Yn[1] += -0.0478571280485400 * F[2][j];
	    Yn[2] += 0.189036518170060 * F[2][j];
	    Yn[3] += 0.328844319980060 * F[2][j];
	    Yn[0] += -0.00990467650730000 * F[3][j];
	    Yn[1] += 0.0160474228065200 * F[3][j];
	    Yn[2] += -0.0241821048998300 * F[3][j];
	    Yn[3] += 0.0625000000000000 * F[3][j];
	    Y[0][j] = Yn[0] * h + y[j];
	    Y[1][j] = Yn[1] * h + y[j];
	    Y[2][j] = Yn[2] * h + y[j];
	    Y[3][j] = Yn[3] * h + y[j];
	  }
#ifdef INSTRUMENT
#pragma omp barrier
	if (me == 0)
	  {
	    double T = time_snap_stop (&time);
#ifdef _OPENMP
	    printf ("#Kernel=23\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=23\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
#pragma omp barrier
    }
//RHS_Approx_Update %8
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
	y[j] =
	  h * 0.220462211176770 *
	  (eval_component (j, t + 0.0885879595126800 * h, Y[0]));
      }
#pragma nounroll_and_jam
    for (int i = 1; i < 4; ++i)
      {
#pragma ivdep
	for (int j = first; j <= last; ++j)
	  {
	    y[j] += h * b[i] * (eval_component (j, t + c[i] * h, Y[i]));
	  }
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