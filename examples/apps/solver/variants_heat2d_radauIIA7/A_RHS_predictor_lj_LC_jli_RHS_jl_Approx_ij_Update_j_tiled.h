#pragma once

#define VARIANT_ID 123456

#include <math.h>
#include "RHS_Heat2D.h"
#include "ODE_radauIIA7.h"
#include "DS_A.h"
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
  for (int k = 0; k < 6; ++k)
    {
#pragma omp barrier
//LC %20
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
	      }
	  }
#ifdef INSTRUMENT
#pragma omp barrier
	if (me == 0)
	  {
	    double T = time_snap_stop (&time);
#ifdef _OPENMP
	    printf ("#Kernel=20\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=20\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
#pragma omp barrier
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
    }
#pragma omp barrier
//Approx %9
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
	 imin (B_Approx_ij, last - first + 1), limit_Approx_ij =
	 imax (first, last + 1 - B_Approx_ij); jj <= last;
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
	printf ("#Kernel=9\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=9\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
//Update %23
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
	 imin (B_Update_j, last - first + 1), limit_Update_j =
	 imax (first, last + 1 - B_Update_j); jj <= last;
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
	printf ("#Kernel=23\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=23\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
}
