#pragma once

#define VARIANT_ID 123456

#include <math.h>
#include "RHS_Cusp.h"
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
      printf ("\n#ImplVariant-123456\n");
    }
#endif
#pragma omp barrier
//RHS_LC_predictor %18
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
	 imin (B_RHS_LC_predictor_jli, last - first + 1),
	 limit_RHS_LC_predictor_jli =
	 imax (first, last + 1 - B_RHS_LC_predictor_jli); jj <= last;
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
	printf ("#Kernel=18\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=18\t#Threads=1\t%.20e\n", T / 1e9 / n);
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

//RHS_LC %31
#ifdef INSTRUMENT
      {
#pragma omp barrier
	time_snap_t time;
	if (me == 0)
	  {
	    time_snap_start (&time);
	  }
#endif
	int bs_RHS_LC_ilj, limit_RHS_LC_ilj;
	for (int jj = first, bs_RHS_LC_ilj =
	     imin (B_RHS_LC_ilj, last - first + 1), limit_RHS_LC_ilj =
	     imax (first, last + 1 - B_RHS_LC_ilj); jj <= last;
	     bs_RHS_LC_ilj = last + 1 - jj, limit_RHS_LC_ilj = last)
	  {
	    for (; jj <= limit_RHS_LC_ilj; jj += bs_RHS_LC_ilj)
	      {
		eval_range (first, last, t + 0.0885879595126800 * h, Yprev[0],
			    Fn);
		for (int l = 0; l < 4; ++l)
		  {
#pragma ivdep
		    for (int j = jj; j < bs_RHS_LC_ilj + jj; ++j)
		      {
			Ycur[l][j] = A[l][0] * Fn[j];
		      }
		  }
#pragma nounroll_and_jam
		for (int i = 1; i < 4; ++i)
		  {
		    eval_range (first, last, t + c[i] * h, Yprev[i], Fn);
		    for (int l = 0; l < 4; ++l)
		      {
#pragma ivdep
			for (int j = jj; j < bs_RHS_LC_ilj + jj; ++j)
			  {
			    Ycur[l][j] += A[l][i] * Fn[j];
			  }
		      }
		  }
#pragma nounroll_and_jam
		for (int l = 0; l < 4; ++l)
		  {
#pragma ivdep
		    for (int j = jj; j < bs_RHS_LC_ilj + jj; ++j)
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
	    printf ("#Kernel=31\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=31\t#Threads=1\t%.20e\n", T / 1e9 / n);
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
	printf ("#Kernel=27\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=27\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
}
