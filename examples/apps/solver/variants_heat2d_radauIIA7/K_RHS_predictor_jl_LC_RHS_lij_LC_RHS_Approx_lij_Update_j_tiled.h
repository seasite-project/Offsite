#pragma once

#define VARIANT_ID 123456

#include <math.h>
#include "RHS_Heat2D.h"
#include "ODE_radauIIA7.h"
#include "DS_K.h"
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
//RHS_predictor %12
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_RHS_predictor_jl, limit_RHS_predictor_jl;
    for (int jj = first, bs_RHS_predictor_jl =
	 imin (B_RHS_predictor_jl, last - first + 1), limit_RHS_predictor_jl =
	 imax (first, last + 1 - B_RHS_predictor_jl); jj <= last;
	 bs_RHS_predictor_jl = last + 1 - jj, limit_RHS_predictor_jl = last)
      {
	for (; jj <= limit_RHS_predictor_jl; jj += bs_RHS_predictor_jl)
	  {
#pragma ivdep
	    for (int j = jj; j < bs_RHS_predictor_jl + jj; ++j)
	      {
		F[0][j] = eval_component (j, t + 0.0885879595126800 * h, y);
		F[1][j] = eval_component (j, t + 0.409466864440740 * h, y);
		F[2][j] = eval_component (j, t + 0.787659461760850 * h, y);
		F[3][j] = eval_component (j, t + 1.00000000000000 * h, y);
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
	double **swp_tmp = Fprev;
	Fprev = F;
	F = swp_tmp;
      }

//LC_RHS %30
#ifdef INSTRUMENT
      {
#pragma omp barrier
	time_snap_t time;
	if (me == 0)
	  {
	    time_snap_start (&time);
	  }
#endif
	int bs_LC_RHS_lij, limit_LC_RHS_lij;
	for (int jj = first, bs_LC_RHS_lij =
	     imin (B_LC_RHS_lij, last - first + 1), limit_LC_RHS_lij =
	     imax (first, last + 1 - B_LC_RHS_lij); jj <= last;
	     bs_LC_RHS_lij = last + 1 - jj, limit_LC_RHS_lij = last)
	  {
	    for (; jj <= limit_LC_RHS_lij; jj += bs_LC_RHS_lij)
	      {
#pragma nounroll_and_jam
		for (int l = 0; l < 4; ++l)
		  {
#pragma ivdep
		    for (int j = jj; j < bs_LC_RHS_lij + jj; ++j)
		      {
			Y[j] = A[l][0] * Fprev[0][j];
		      }
		    for (int i = 1; i < 4; ++i)
		      {
#pragma ivdep
			for (int j = jj; j < bs_LC_RHS_lij + jj; ++j)
			  {
			    Y[j] += A[l][i] * Fprev[i][j];
			  }
		      }
#pragma ivdep
		    for (int j = jj; j < bs_LC_RHS_lij + jj; ++j)
		      {
			Y[j] = Y[j] * h + y[j];
		      }
#pragma omp barrier
		    eval_range (first, last, t + c[l] * h, Y, F[l]);
#pragma omp barrier
		  }
	      }
	  }
#ifdef INSTRUMENT
#pragma omp barrier
	if (me == 0)
	  {
	    double T = time_snap_stop (&time);
#ifdef _OPENMP
	    printf ("#Kernel=30\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=30\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
    }

#pragma omp master
  {
    double **swp_tmp = Fprev;
    Fprev = F;
    F = swp_tmp;
  }

//LC_RHS_Approx %3
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_LC_RHS_Approx_lij, limit_LC_RHS_Approx_lij;
    for (int jj = first, bs_LC_RHS_Approx_lij =
	 imin (B_LC_RHS_Approx_lij, last - first + 1),
	 limit_LC_RHS_Approx_lij =
	 imax (first, last + 1 - B_LC_RHS_Approx_lij); jj <= last;
	 bs_LC_RHS_Approx_lij = last + 1 - jj, limit_LC_RHS_Approx_lij = last)
      {
	for (; jj <= limit_LC_RHS_Approx_lij; jj += bs_LC_RHS_Approx_lij)
	  {
#pragma nounroll_and_jam
	    for (int l = 0; l < 4; ++l)
	      {
#pragma ivdep
		for (int j = jj; j < bs_LC_RHS_Approx_lij + jj; ++j)
		  {
		    Y[j] = A[l][0] * Fprev[0][j];
		  }
		for (int i = 1; i < 4; ++i)
		  {
#pragma ivdep
		    for (int j = jj; j < bs_LC_RHS_Approx_lij + jj; ++j)
		      {
			Y[j] += A[l][i] * Fprev[i][j];
		      }
		  }
#pragma ivdep
		for (int j = jj; j < bs_LC_RHS_Approx_lij + jj; ++j)
		  {
		    Y[j] = Y[j] * h + y[j];
		  }
#pragma omp barrier
#pragma ivdep
		for (int j = jj; j < bs_LC_RHS_Approx_lij + jj; ++j)
		  {
		    dy[j] += b[l] * (eval_component (j, t + c[l] * h, Y));
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
	printf ("#Kernel=3\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=3\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
#pragma omp barrier
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
