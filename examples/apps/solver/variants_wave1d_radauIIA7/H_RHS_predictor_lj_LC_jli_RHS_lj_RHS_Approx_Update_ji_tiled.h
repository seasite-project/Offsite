#pragma once

#define VARIANT_ID 123456

#include <math.h>
#include "RHS_Wave1D.h"
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
//RHS_predictor %15
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
	printf ("#Kernel=15\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=15\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
#pragma omp barrier
//LC %24
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_LC_jli, limit_LC_jli;
    for (int jj = first, bs_LC_jli =
	 imin (B_LC_jli, last - first + 1), limit_LC_jli =
	 imax (first, last + 1 - B_LC_jli); jj <= last;
	 bs_LC_jli = last + 1 - jj, limit_LC_jli = last)
      {
	for (; jj <= limit_LC_jli; jj += bs_LC_jli)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_LC_jli + jj; ++j)
	      {
		Y[0][j] = 0.112999479323160 * F[0][j];
		Y[0][j] += -0.0403092207235200 * F[1][j];
		Y[0][j] += 0.0258023774203400 * F[2][j];
		Y[0][j] += -0.00990467650730000 * F[3][j];
		Y[0][j] = Y[0][j] * h + y[j];
		Y[1][j] = 0.234383995747400 * F[0][j];
		Y[1][j] += 0.206892573935360 * F[1][j];
		Y[1][j] += -0.0478571280485400 * F[2][j];
		Y[1][j] += 0.0160474228065200 * F[3][j];
		Y[1][j] = Y[1][j] * h + y[j];
		Y[2][j] = 0.216681784623250 * F[0][j];
		Y[2][j] += 0.406123263867370 * F[1][j];
		Y[2][j] += 0.189036518170060 * F[2][j];
		Y[2][j] += -0.0241821048998300 * F[3][j];
		Y[2][j] = Y[2][j] * h + y[j];
		Y[3][j] = 0.220462211176770 * F[0][j];
		Y[3][j] += 0.388193468843170 * F[1][j];
		Y[3][j] += 0.328844319980060 * F[2][j];
		Y[3][j] += 0.0625000000000000 * F[3][j];
		Y[3][j] = Y[3][j] * h + y[j];
	      }
	  }
      }
#ifdef INSTRUMENT
#pragma omp barrier
    if (me == 0)
      {
	double T = time_snap_stop (&time);
#ifdef _OPENMP
	printf ("#Kernel=24\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=24\t#Threads=1\t%.20e\n", T / 1e9 / n);
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
	int bs_RHS_lj, limit_RHS_lj;
	for (int jj = first, bs_RHS_lj =
	     imin (B_RHS_lj, last - first + 1), limit_RHS_lj =
	     imax (first, last + 1 - B_RHS_lj); jj <= last;
	     bs_RHS_lj = last + 1 - jj, limit_RHS_lj = last)
	  {
	    for (; jj <= limit_RHS_lj; jj += bs_RHS_lj)
	      {
		for (int l = 0; l < 4; ++l)
		  {
		    eval_range (first, last, t + c[l] * h, Y[l], F[l]);
		  }
	      }
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
//LC %24
#ifdef INSTRUMENT
      {
#pragma omp barrier
	time_snap_t time;
	if (me == 0)
	  {
	    time_snap_start (&time);
	  }
#endif
	int bs_LC_jli, limit_LC_jli;
	for (int jj = first, bs_LC_jli =
	     imin (B_LC_jli, last - first + 1), limit_LC_jli =
	     imax (first, last + 1 - B_LC_jli); jj <= last;
	     bs_LC_jli = last + 1 - jj, limit_LC_jli = last)
	  {
	    for (; jj <= limit_LC_jli; jj += bs_LC_jli)
	      {
#pragma ivdep
		for (int j = jj; j < bs_LC_jli + jj; ++j)
		  {
		    Y[0][j] = 0.112999479323160 * F[0][j];
		    Y[0][j] += -0.0403092207235200 * F[1][j];
		    Y[0][j] += 0.0258023774203400 * F[2][j];
		    Y[0][j] += -0.00990467650730000 * F[3][j];
		    Y[0][j] = Y[0][j] * h + y[j];
		    Y[1][j] = 0.234383995747400 * F[0][j];
		    Y[1][j] += 0.206892573935360 * F[1][j];
		    Y[1][j] += -0.0478571280485400 * F[2][j];
		    Y[1][j] += 0.0160474228065200 * F[3][j];
		    Y[1][j] = Y[1][j] * h + y[j];
		    Y[2][j] = 0.216681784623250 * F[0][j];
		    Y[2][j] += 0.406123263867370 * F[1][j];
		    Y[2][j] += 0.189036518170060 * F[2][j];
		    Y[2][j] += -0.0241821048998300 * F[3][j];
		    Y[2][j] = Y[2][j] * h + y[j];
		    Y[3][j] = 0.220462211176770 * F[0][j];
		    Y[3][j] += 0.388193468843170 * F[1][j];
		    Y[3][j] += 0.328844319980060 * F[2][j];
		    Y[3][j] += 0.0625000000000000 * F[3][j];
		    Y[3][j] = Y[3][j] * h + y[j];
		  }
	      }
	  }
#ifdef INSTRUMENT
#pragma omp barrier
	if (me == 0)
	  {
	    double T = time_snap_stop (&time);
#ifdef _OPENMP
	    printf ("#Kernel=24\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=24\t#Threads=1\t%.20e\n", T / 1e9 / n);
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
		y[j] += h * dy[j];
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
