#pragma once

#define VARIANT_ID 195

#include <math.h>
#include "RHS_InverterChain.h"
#include "ODE_radauIIA7.h"
#include "DS_C.h"
#ifdef INSTRUMENT
#include "timesnap.h"
#endif

void
timestep (const int me, const int first, const int last, double t, double h)
{
#ifdef INSTRUMENT
  if (me == 0)
    {
      printf ("\n#ImplVariant-195\n");
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
    int bs_RHS_LC_predictor_ilj, limit_RHS_LC_predictor_ilj;
    for (int jj = first, bs_RHS_LC_predictor_ilj =
	 imin (B, last - first + 1), limit_RHS_LC_predictor_ilj =
	 imax (first, last + 1 - B); jj <= last;
	 bs_RHS_LC_predictor_ilj = last + 1 - jj, limit_RHS_LC_predictor_ilj =
	 last)
      {
	for (; jj <= limit_RHS_LC_predictor_ilj;
	     jj += bs_RHS_LC_predictor_ilj)
	  {
	    eval_range (first, last, t + 0.0885879595126800 * h, y, Fn);
	    for (int l = 0; l < 4; ++l)
	      {
#pragma ivdep
		for (int j = jj; j < bs_RHS_LC_predictor_ilj + jj; ++j)
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
		    for (int j = jj; j < bs_RHS_LC_predictor_ilj + jj; ++j)
		      {
			Ycur[l][j] += A[l][i] * Fn[j];
		      }
		  }
	      }
#pragma nounroll_and_jam
	    for (int l = 0; l < 4; ++l)
	      {
#pragma ivdep
		for (int j = jj; j < bs_RHS_LC_predictor_ilj + jj; ++j)
		  {
		    Ycur[l][j] = Ycur[l][j] * h + y[j];
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
	     imin (B, last - first + 1), limit_RHS_LC_jil =
	     imax (first, last + 1 - B); jj <= last;
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
		    tmp[0] = 0.112999479323160 * f;
		    tmp[1] = 0.234383995747400 * f;
		    tmp[2] = 0.216681784623250 * f;
		    tmp[3] = 0.220462211176770 * f;
		    f =
		      eval_component (j, t + 0.409466864440740 * h, Yprev[1]);
		    tmp[0] += -0.0403092207235200 * f;
		    tmp[1] += 0.206892573935360 * f;
		    tmp[2] += 0.406123263867370 * f;
		    tmp[3] += 0.388193468843170 * f;
		    f =
		      eval_component (j, t + 0.787659461760850 * h, Yprev[2]);
		    tmp[0] += 0.0258023774203400 * f;
		    tmp[1] += -0.0478571280485400 * f;
		    tmp[2] += 0.189036518170060 * f;
		    tmp[3] += 0.328844319980060 * f;
		    f =
		      eval_component (j, t + 1.00000000000000 * h, Yprev[3]);
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

//RHS %33
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
	 imin (B, last - first + 1), limit_RHS_jl =
	 imax (first, last + 1 - B); jj <= last;
	 bs_RHS_jl = last + 1 - jj, limit_RHS_jl = last)
      {
	for (; jj <= limit_RHS_jl; jj += bs_RHS_jl)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_RHS_jl + jj; ++j)
	      {
		F[0][j] =
		  eval_component (j, t + 0.0885879595126800 * h, Yprev[0]);
		F[1][j] =
		  eval_component (j, t + 0.409466864440740 * h, Yprev[1]);
		F[2][j] =
		  eval_component (j, t + 0.787659461760850 * h, Yprev[2]);
		F[3][j] =
		  eval_component (j, t + 1.00000000000000 * h, Yprev[3]);
	      }
	  }
      }
#ifdef INSTRUMENT
#pragma omp barrier
    if (me == 0)
      {
	double T = time_snap_stop (&time);
#ifdef _OPENMP
	printf ("#Kernel=33\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=33\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
#pragma omp barrier
//Approx %13
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_Approx_ij, limit_Approx_ij;
    for (int jj = first, bs_Approx_ij =
	 imin (B, last - first + 1), limit_Approx_ij =
	 imax (first, last + 1 - B); jj <= last;
	 bs_Approx_ij = last + 1 - jj, limit_Approx_ij = last)
      {
	for (; jj <= limit_Approx_ij; jj += bs_Approx_ij)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_Approx_ij + jj; ++j)
	      {
		dy[j] = 0.220462211176770 * F[0][j];
	      }
#pragma nounroll_and_jam
	    for (int i = 1; i < 4; ++i)
	      {
#pragma ivdep
		for (int j = jj; j < bs_Approx_ij + jj; ++j)
		  {
		    dy[j] += b[i] * F[i][j];
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
	printf ("#Kernel=13\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=13\t#Threads=1\t%.20e\n", T / 1e9 / n);
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
	printf ("#Kernel=27\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=27\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
}
