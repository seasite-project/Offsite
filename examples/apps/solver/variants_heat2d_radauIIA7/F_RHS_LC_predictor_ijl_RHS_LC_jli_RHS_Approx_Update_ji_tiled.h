#pragma once

#define VARIANT_ID 123456

#include <math.h>
#include "RHS_Heat2D.h"
#include "ODE_radauIIA7.h"
#include "DS_F.h"
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
    int bs_RHS_LC_predictor_ijl, limit_RHS_LC_predictor_ijl;
    for (int jj = first, bs_RHS_LC_predictor_ijl =
	 imin (B_RHS_LC_predictor_ijl, last - first + 1),
	 limit_RHS_LC_predictor_ijl =
	 imax (first, last + 1 - B_RHS_LC_predictor_ijl); jj <= last;
	 bs_RHS_LC_predictor_ijl = last + 1 - jj, limit_RHS_LC_predictor_ijl =
	 last)
      {
	for (; jj <= limit_RHS_LC_predictor_ijl;
	     jj += bs_RHS_LC_predictor_ijl)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_RHS_LC_predictor_ijl + jj; ++j)
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
		for (int j = jj; j < bs_RHS_LC_predictor_ijl + jj; ++j)
		  {
		    f = eval_component (j, t + c[i] * h, y);
		    Ycur[0][j] += A[0][i] * f;
		    Ycur[1][j] += A[1][i] * f;
		    Ycur[2][j] += A[2][i] * f;
		    Ycur[3][j] += A[3][i] * f;
		  }
	      }
#pragma ivdep
	    for (int j = jj; j < bs_RHS_LC_predictor_ijl + jj; ++j)
	      {
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
	int bs_RHS_LC_jli, limit_RHS_LC_jli;
	for (int jj = first, bs_RHS_LC_jli =
	     imin (B_RHS_LC_jli, last - first + 1), limit_RHS_LC_jli =
	     imax (first, last + 1 - B_RHS_LC_jli); jj <= last;
	     bs_RHS_LC_jli = last + 1 - jj, limit_RHS_LC_jli = last)
	  {
	    for (; jj <= limit_RHS_LC_jli; jj += bs_RHS_LC_jli)
	      {
#pragma ivdep
		for (int j = jj; j < bs_RHS_LC_jli + jj; ++j)
		  {
		    Fs[0] =
		      eval_component (j, t + 0.0885879595126800 * h,
				      Yprev[0]);
		    Fs[1] =
		      eval_component (j, t + 0.409466864440740 * h, Yprev[1]);
		    Fs[2] =
		      eval_component (j, t + 0.787659461760850 * h, Yprev[2]);
		    Fs[3] =
		      eval_component (j, t + 1.00000000000000 * h, Yprev[3]);
		    Ycur[0][j] = 0.112999479323160 * Fs[0];
		    Ycur[0][j] += -0.0403092207235200 * Fs[0];
		    Ycur[0][j] += 0.0258023774203400 * Fs[0];
		    Ycur[0][j] += -0.00990467650730000 * Fs[0];
		    Ycur[0][j] = Ycur[0][j] * h + y[j];
		    Ycur[1][j] = 0.234383995747400 * Fs[1];
		    Ycur[1][j] += 0.206892573935360 * Fs[1];
		    Ycur[1][j] += -0.0478571280485400 * Fs[1];
		    Ycur[1][j] += 0.0160474228065200 * Fs[1];
		    Ycur[1][j] = Ycur[1][j] * h + y[j];
		    Ycur[2][j] = 0.216681784623250 * Fs[2];
		    Ycur[2][j] += 0.406123263867370 * Fs[2];
		    Ycur[2][j] += 0.189036518170060 * Fs[2];
		    Ycur[2][j] += -0.0241821048998300 * Fs[2];
		    Ycur[2][j] = Ycur[2][j] * h + y[j];
		    Ycur[3][j] = 0.220462211176770 * Fs[3];
		    Ycur[3][j] += 0.388193468843170 * Fs[3];
		    Ycur[3][j] += 0.328844319980060 * Fs[3];
		    Ycur[3][j] += 0.0625000000000000 * Fs[3];
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
    double **tmp = Yprev;
    Yprev = Ycur;
    Ycur = tmp;
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
