#pragma once

#define VARIANT_ID 123456

#include <math.h>
#include "RHS_Cusp.h"
#include "ODE_radauIIA7.h"
#include "DS_D.h"
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
//RHS_LC_predictor %17
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_RHS_LC_predictor_jil, limit_RHS_LC_predictor_jil;
    for (int jj = first, bs_RHS_LC_predictor_jil =
	 imin (B_RHS_LC_predictor_jil, last - first + 1),
	 limit_RHS_LC_predictor_jil =
	 imax (first, last + 1 - B_RHS_LC_predictor_jil); jj <= last;
	 bs_RHS_LC_predictor_jil = last + 1 - jj, limit_RHS_LC_predictor_jil =
	 last)
      {
	for (; jj <= limit_RHS_LC_predictor_jil;
	     jj += bs_RHS_LC_predictor_jil)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_RHS_LC_predictor_jil + jj; ++j)
	      {
		f = eval_component (j, t + 0.0885879595126800 * h, y);
		Ycur[0][j] = 0.112999479323160 * f;
		Ycur[1][j] = 0.234383995747400 * f;
		Ycur[2][j] = 0.216681784623250 * f;
		Ycur[3][j] = 0.220462211176770 * f;
		f = eval_component (j, t + 0.409466864440740 * h, y);
		Ycur[0][j] += -0.0403092207235200 * f;
		Ycur[1][j] += 0.206892573935360 * f;
		Ycur[2][j] += 0.406123263867370 * f;
		Ycur[3][j] += 0.388193468843170 * f;
		f = eval_component (j, t + 0.787659461760850 * h, y);
		Ycur[0][j] += 0.0258023774203400 * f;
		Ycur[1][j] += -0.0478571280485400 * f;
		Ycur[2][j] += 0.189036518170060 * f;
		Ycur[3][j] += 0.328844319980060 * f;
		f = eval_component (j, t + 1.00000000000000 * h, y);
		Ycur[0][j] += -0.00990467650730000 * f;
		Ycur[1][j] += 0.0160474228065200 * f;
		Ycur[2][j] += -0.0241821048998300 * f;
		Ycur[3][j] += 0.0625000000000000 * f;
		Ycur[0][j] = Ycur[0][j] * h + y[j];
		Ycur[1][j] = Ycur[1][j] * h + y[j];
		Ycur[2][j] = Ycur[2][j] * h + y[j];
		Ycur[3][j] = Ycur[3][j] * h + y[j];
	      }
	  }
      }
#ifdef INSTRUMENT
#pragma omp barrier
    if (me == 0)
      {
	double T = time_snap_stop (&time);
#ifdef _OPENMP
	printf ("#Kernel=17\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=17\t#Threads=1\t%.20e\n", T / 1e9 / n);
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
	int bs_RHS_LC_jil, limit_RHS_LC_jil;
	for (int jj = first, bs_RHS_LC_jil =
	     imin (B_RHS_LC_jil, last - first + 1), limit_RHS_LC_jil =
	     imax (first, last + 1 - B_RHS_LC_jil); jj <= last;
	     bs_RHS_LC_jil = last + 1 - jj, limit_RHS_LC_jil = last)
	  {
	    for (; jj <= limit_RHS_LC_jil; jj += bs_RHS_LC_jil)
	      {
#pragma ivdep
		for (int j = jj; j < bs_RHS_LC_jil + jj; ++j)
		  {
		    f =
		      eval_component (j, t + 0.0885879595126800 * h,
				      Yprev[0]);
		    Ycur[0][j] = 0.112999479323160 * f;
		    Ycur[1][j] = 0.234383995747400 * f;
		    Ycur[2][j] = 0.216681784623250 * f;
		    Ycur[3][j] = 0.220462211176770 * f;
		    f =
		      eval_component (j, t + 0.409466864440740 * h, Yprev[1]);
		    Ycur[0][j] += -0.0403092207235200 * f;
		    Ycur[1][j] += 0.206892573935360 * f;
		    Ycur[2][j] += 0.406123263867370 * f;
		    Ycur[3][j] += 0.388193468843170 * f;
		    f =
		      eval_component (j, t + 0.787659461760850 * h, Yprev[2]);
		    Ycur[0][j] += 0.0258023774203400 * f;
		    Ycur[1][j] += -0.0478571280485400 * f;
		    Ycur[2][j] += 0.189036518170060 * f;
		    Ycur[3][j] += 0.328844319980060 * f;
		    f =
		      eval_component (j, t + 1.00000000000000 * h, Yprev[3]);
		    Ycur[0][j] += -0.00990467650730000 * f;
		    Ycur[1][j] += 0.0160474228065200 * f;
		    Ycur[2][j] += -0.0241821048998300 * f;
		    Ycur[3][j] += 0.0625000000000000 * f;
		    Ycur[0][j] = Ycur[0][j] * h + y[j];
		    Ycur[1][j] = Ycur[1][j] * h + y[j];
		    Ycur[2][j] = Ycur[2][j] * h + y[j];
		    Ycur[3][j] = Ycur[3][j] * h + y[j];
		  }
	      }
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
		eval_range (first, last, t + c[l] * h, Yprev[l], F[l]);
	      }
	  }
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
//ApproxUpdate %5
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_ApproxUpdate_ij, limit_ApproxUpdate_ij;
    for (int jj = first, bs_ApproxUpdate_ij =
	 imin (B_ApproxUpdate_ij, last - first + 1), limit_ApproxUpdate_ij =
	 imax (first, last + 1 - B_ApproxUpdate_ij); jj <= last;
	 bs_ApproxUpdate_ij = last + 1 - jj, limit_ApproxUpdate_ij = last)
      {
	for (; jj <= limit_ApproxUpdate_ij; jj += bs_ApproxUpdate_ij)
	  {
#pragma nounroll_and_jam
	    for (int i = 0; i < 4; ++i)
	      {
#pragma ivdep
		for (int j = jj; j < bs_ApproxUpdate_ij + jj; ++j)
		  {
		    y[j] += h * b[i] * F[i][j];
		  }
	      }
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
}
