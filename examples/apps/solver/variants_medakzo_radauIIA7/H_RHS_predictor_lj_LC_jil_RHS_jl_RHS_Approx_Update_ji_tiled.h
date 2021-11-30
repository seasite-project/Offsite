#pragma once

#define VARIANT_ID 123456

#include <math.h>
#include "RHS_Medakzo.h"
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
      printf ("\n#ImplVariant-123456\n");
    }
#endif
#pragma omp barrier
//RHS_predictor %11
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_RHS_predictor_lj, limit_RHS_predictor_lj;
    for (int jj = first, bs_RHS_predictor_lj =
	 imin (B_RHS_predictor_lj, last - first + 1), limit_RHS_predictor_lj =
	 imax (first, last + 1 - B_RHS_predictor_lj); jj <= last;
	 bs_RHS_predictor_lj = last + 1 - jj, limit_RHS_predictor_lj = last)
      {
	for (; jj <= limit_RHS_predictor_lj; jj += bs_RHS_predictor_lj)
	  {
	    for (int l = 0; l < 4; ++l)
	      {
		eval_range (first, last, t + c[l] * h, y, F[l]);
	      }
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
#pragma omp barrier
//LC %19
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_LC_jil, limit_LC_jil;
    for (int jj = first, bs_LC_jil =
	 imin (B_LC_jil, last - first + 1), limit_LC_jil =
	 imax (first, last + 1 - B_LC_jil); jj <= last;
	 bs_LC_jil = last + 1 - jj, limit_LC_jil = last)
      {
	for (; jj <= limit_LC_jil; jj += bs_LC_jil)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_LC_jil + jj; ++j)
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
	  }
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
//RHS %29
#ifdef INSTRUMENT
      {
#pragma omp barrier
	time_snap_t time;
	if (me == 0)
	  {
	    time_snap_start (&time);
	  }
#endif
	int bs_RHS_jl, limit_RHS_jl;
	for (int jj = first, bs_RHS_jl =
	     imin (B_RHS_jl, last - first + 1), limit_RHS_jl =
	     imax (first, last + 1 - B_RHS_jl); jj <= last;
	     bs_RHS_jl = last + 1 - jj, limit_RHS_jl = last)
	  {
	    for (; jj <= limit_RHS_jl; jj += bs_RHS_jl)
	      {
#pragma ivdep
		for (int j = jj; j < bs_RHS_jl + jj; ++j)
		  {
		    F[0][j] =
		      eval_component (j, t + 0.0885879595126800 * h, Y[0]);
		    F[1][j] =
		      eval_component (j, t + 0.409466864440740 * h, Y[1]);
		    F[2][j] =
		      eval_component (j, t + 0.787659461760850 * h, Y[2]);
		    F[3][j] =
		      eval_component (j, t + 1.00000000000000 * h, Y[3]);
		  }
	      }
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
//LC %19
#ifdef INSTRUMENT
      {
#pragma omp barrier
	time_snap_t time;
	if (me == 0)
	  {
	    time_snap_start (&time);
	  }
#endif
	int bs_LC_jil, limit_LC_jil;
	for (int jj = first, bs_LC_jil =
	     imin (B_LC_jil, last - first + 1), limit_LC_jil =
	     imax (first, last + 1 - B_LC_jil); jj <= last;
	     bs_LC_jil = last + 1 - jj, limit_LC_jil = last)
	  {
	    for (; jj <= limit_LC_jil; jj += bs_LC_jil)
	      {
#pragma ivdep
		for (int j = jj; j < bs_LC_jil + jj; ++j)
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
	      }
	  }
#ifdef INSTRUMENT
#pragma omp barrier
	if (me == 0)
	  {
	    double T = time_snap_stop (&time);
#ifdef _OPENMP
	    printf ("#Kernel=19\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=19\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
#pragma omp barrier
    }
//RHS_Approx_Update %1
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_RHS_Approx_Update_ji, limit_RHS_Approx_Update_ji;
    for (int jj = first, bs_RHS_Approx_Update_ji =
	 imin (B_RHS_Approx_Update_ji, last - first + 1),
	 limit_RHS_Approx_Update_ji =
	 imax (first, last + 1 - B_RHS_Approx_Update_ji); jj <= last;
	 bs_RHS_Approx_Update_ji = last + 1 - jj, limit_RHS_Approx_Update_ji =
	 last)
      {
	for (; jj <= limit_RHS_Approx_Update_ji;
	     jj += bs_RHS_Approx_Update_ji)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_RHS_Approx_Update_ji + jj; ++j)
	      {
		dy =
		  0.220462211176770 *
		  (eval_component (j, t + 0.0885879595126800 * h, Y[0]));
		dy +=
		  0.388193468843170 *
		  (eval_component (j, t + 0.409466864440740 * h, Y[1]));
		dy +=
		  0.328844319980060 *
		  (eval_component (j, t + 0.787659461760850 * h, Y[2]));
		dy +=
		  0.0625000000000000 *
		  (eval_component (j, t + 1.00000000000000 * h, Y[3]));
		y[j] += h * dy;
	      }
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
}
