#pragma once

#define VARIANT_ID 332

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
      printf ("\n#ImplVariant-332\n");
    }
#endif
#pragma omp barrier
//RHS_LC_predictor %20
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
	printf ("#Kernel=20\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=20\t#Threads=1\t%.20e\n", T / 1e9 / n);
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

//RHS_LC %29
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
	    Fs[0] = eval_component (j, t + 0.0885879595126800 * h, Yprev[0]);
	    Fs[1] = eval_component (j, t + 0.409466864440740 * h, Yprev[1]);
	    Fs[2] = eval_component (j, t + 0.787659461760850 * h, Yprev[2]);
	    Fs[3] = eval_component (j, t + 1.00000000000000 * h, Yprev[3]);
	    tmp_s = 0.112999479323160 * Fs[0];
	    tmp_s += -0.0403092207235200 * Fs[0];
	    tmp_s += 0.0258023774203400 * Fs[0];
	    tmp_s += -0.00990467650730000 * Fs[0];
	    Ycur[0][j] = tmp_s * h + y[j];
	    tmp_s = 0.234383995747400 * Fs[1];
	    tmp_s += 0.206892573935360 * Fs[1];
	    tmp_s += -0.0478571280485400 * Fs[1];
	    tmp_s += 0.0160474228065200 * Fs[1];
	    Ycur[1][j] = tmp_s * h + y[j];
	    tmp_s = 0.216681784623250 * Fs[2];
	    tmp_s += 0.406123263867370 * Fs[2];
	    tmp_s += 0.189036518170060 * Fs[2];
	    tmp_s += -0.0241821048998300 * Fs[2];
	    Ycur[2][j] = tmp_s * h + y[j];
	    tmp_s = 0.220462211176770 * Fs[3];
	    tmp_s += 0.388193468843170 * Fs[3];
	    tmp_s += 0.328844319980060 * Fs[3];
	    tmp_s += 0.0625000000000000 * Fs[3];
	    Ycur[3][j] = tmp_s * h + y[j];
	  }
#ifdef INSTRUMENT
#pragma omp barrier
	if (me == 0)
	  {
	    double T = time_snap_stop (&time);
#ifdef _OPENMP
	    printf ("#Kernel=29\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=29\t#Threads=1\t%.20e\n", T / 1e9 / n);
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

//RHS_Approx %12
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
	dy[j] +=
	  0.388193468843170 *
	  (eval_component (j, t + 0.409466864440740 * h, Yprev[1]));
	dy[j] +=
	  0.328844319980060 *
	  (eval_component (j, t + 0.787659461760850 * h, Yprev[2]));
	dy[j] +=
	  0.0625000000000000 *
	  (eval_component (j, t + 1.00000000000000 * h, Yprev[3]));
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
