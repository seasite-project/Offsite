#pragma once

#define VARIANT_ID 307

#include <math.h>
#include "RHS_Heat2D.h"
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
      printf ("\n#ImplVariant-307\n");
    }
#endif
#pragma omp barrier
//RHS_LC_predictor %11
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
	 imin (B, last - first + 1), limit_RHS_LC_predictor_jil =
	 imax (first, last + 1 - B); jj <= last;
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
		tmp[0] = 0.112999479323160 * f;
		tmp[1] = 0.234383995747400 * f;
		tmp[2] = 0.216681784623250 * f;
		tmp[3] = 0.220462211176770 * f;
		f = eval_component (j, t + 0.409466864440740 * h, y);
		tmp[0] += -0.0403092207235200 * f;
		tmp[1] += 0.206892573935360 * f;
		tmp[2] += 0.406123263867370 * f;
		tmp[3] += 0.388193468843170 * f;
		f = eval_component (j, t + 0.787659461760850 * h, y);
		tmp[0] += 0.0258023774203400 * f;
		tmp[1] += -0.0478571280485400 * f;
		tmp[2] += 0.189036518170060 * f;
		tmp[3] += 0.328844319980060 * f;
		f = eval_component (j, t + 1.00000000000000 * h, y);
		tmp[0] += -0.00990467650730000 * f;
		tmp[1] += 0.0160474228065200 * f;
		tmp[2] += -0.0241821048998300 * f;
		tmp[3] += 0.0625000000000000 * f;
		Ycur[0][j] = tmp[0] * h + y[j];
		Ycur[1][j] = tmp[1] * h + y[j];
		Ycur[2][j] = tmp[2] * h + y[j];
		Ycur[3][j] = tmp[3] * h + y[j];
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
  for (int k = 0; k < 5; ++k)
    {

#pragma omp master
      {
	double **swp_tmp = Yprev;
	Yprev = Ycur;
	Ycur = swp_tmp;
      }

//RHS_LC %23
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
	     imin (B, last - first + 1), limit_RHS_LC_jli =
	     imax (first, last + 1 - B); jj <= last;
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

#pragma omp master
  {
    double **swp_tmp = Yprev;
    Yprev = Ycur;
    Ycur = swp_tmp;
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
    int bs_RHS_Approx_ij, limit_RHS_Approx_ij;
    for (int jj = first, bs_RHS_Approx_ij =
	 imin (B, last - first + 1), limit_RHS_Approx_ij =
	 imax (first, last + 1 - B); jj <= last;
	 bs_RHS_Approx_ij = last + 1 - jj, limit_RHS_Approx_ij = last)
      {
	for (; jj <= limit_RHS_Approx_ij; jj += bs_RHS_Approx_ij)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_RHS_Approx_ij + jj; ++j)
	      {
		dy[j] =
		  0.220462211176770 *
		  (eval_component (j, t + 0.0885879595126800 * h, Yprev[0]));
	      }
	    for (int i = 1; i < 4; ++i)
	      {
#pragma ivdep
		for (int j = jj; j < bs_RHS_Approx_ij + jj; ++j)
		  {
		    dy[j] +=
		      b[i] * (eval_component (j, t + c[i] * h, Yprev[i]));
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
    int bs_Update_j, limit_Update_j;
    for (int jj = first, bs_Update_j =
	 imin (B, last - first + 1), limit_Update_j =
	 imax (first, last + 1 - B); jj <= last;
	 bs_Update_j = last + 1 - jj, limit_Update_j = last)
      {
	for (; jj <= limit_Update_j; jj += bs_Update_j)
	  {
#pragma nounroll_and_jam
#pragma ivdep
	    for (int j = jj; j < bs_Update_j + jj; ++j)
	      {
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
	printf ("#Kernel=21\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=21\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
}
