#pragma once

#define VARIANT_ID 398

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
      printf ("\n#ImplVariant-398\n");
    }
#endif
#pragma omp barrier
//RHS_LC_predictor %12
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_RHS_LC_predictor_jli, limit_RHS_LC_predictor_jli;
    for (int jj = first, bs_RHS_LC_predictor_jli =
	 imin (B, last - first + 1), limit_RHS_LC_predictor_jli =
	 imax (first, last + 1 - B); jj <= last;
	 bs_RHS_LC_predictor_jli = last + 1 - jj, limit_RHS_LC_predictor_jli =
	 last)
      {
	for (; jj <= limit_RHS_LC_predictor_jli;
	     jj += bs_RHS_LC_predictor_jli)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_RHS_LC_predictor_jli + jj; ++j)
	      {
		Fs[0] = eval_component (j, t + 0.0885879595126800 * h, y);
		Fs[1] = eval_component (j, t + 0.409466864440740 * h, y);
		Fs[2] = eval_component (j, t + 0.787659461760850 * h, y);
		Fs[3] = eval_component (j, t + 1.00000000000000 * h, y);
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
	  }
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
#pragma omp barrier
  for (int k = 0; k < 5; ++k)
    {

#pragma omp master
      {
	double **swp_tmp = Yprev;
	Yprev = Ycur;
	Ycur = swp_tmp;
      }

//RHS_LC %24
#ifdef INSTRUMENT
      {
#pragma omp barrier
	time_snap_t time;
	if (me == 0)
	  {
	    time_snap_start (&time);
	  }
#endif
	int bs_RHS_LC_ijl, limit_RHS_LC_ijl;
	for (int jj = first, bs_RHS_LC_ijl =
	     imin (B, last - first + 1), limit_RHS_LC_ijl =
	     imax (first, last + 1 - B); jj <= last;
	     bs_RHS_LC_ijl = last + 1 - jj, limit_RHS_LC_ijl = last)
	  {
	    for (; jj <= limit_RHS_LC_ijl; jj += bs_RHS_LC_ijl)
	      {
#pragma ivdep
		for (int j = jj; j < bs_RHS_LC_ijl + jj; ++j)
		  {
		    f =
		      eval_component (j, t + 0.0885879595126800 * h,
				      Yprev[0]);
		    Ycur[0][j] = 0.112999479323160 * f;
		    Ycur[1][j] = 0.234383995747400 * f;
		    Ycur[2][j] = 0.216681784623250 * f;
		    Ycur[3][j] = 0.220462211176770 * f;
		  }
#pragma nounroll_and_jam
		for (int i = 1; i < 4; ++i)
		  {
#pragma ivdep
		    for (int j = jj; j < bs_RHS_LC_ijl + jj; ++j)
		      {
			f = eval_component (j, t + c[i] * h, Yprev[i]);
			Ycur[0][j] += A[0][i] * f;
			Ycur[1][j] += A[1][i] * f;
			Ycur[2][j] += A[2][i] * f;
			Ycur[3][j] += A[3][i] * f;
		      }
		  }
#pragma ivdep
		for (int j = jj; j < bs_RHS_LC_ijl + jj; ++j)
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

#pragma omp master
  {
    double **swp_tmp = Yprev;
    Yprev = Ycur;
    Ycur = swp_tmp;
  }

//RHS_Approx_Update %2
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_RHS_Approx_Update_ij, limit_RHS_Approx_Update_ij;
    for (int jj = first, bs_RHS_Approx_Update_ij =
	 imin (B, last - first + 1), limit_RHS_Approx_Update_ij =
	 imax (first, last + 1 - B); jj <= last;
	 bs_RHS_Approx_Update_ij = last + 1 - jj, limit_RHS_Approx_Update_ij =
	 last)
      {
	for (; jj <= limit_RHS_Approx_Update_ij;
	     jj += bs_RHS_Approx_Update_ij)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_RHS_Approx_Update_ij + jj; ++j)
	      {
		y[j] =
		  h * 0.220462211176770 *
		  (eval_component (j, t + 0.0885879595126800 * h, Yprev[0]));
	      }
#pragma nounroll_and_jam
	    for (int i = 1; i < 4; ++i)
	      {
#pragma ivdep
		for (int j = jj; j < bs_RHS_Approx_Update_ij + jj; ++j)
		  {
		    y[j] +=
		      h * b[i] * (eval_component (j, t + c[i] * h, Yprev[i]));
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
	printf ("#Kernel=2\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=2\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
}
