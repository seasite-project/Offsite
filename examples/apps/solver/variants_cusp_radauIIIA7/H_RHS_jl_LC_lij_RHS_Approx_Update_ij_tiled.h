#pragma once

#define VARIANT_ID 123456

#include <math.h>
#include "RHS_Cusp.h"
#include "ODE_radauIIA7.h"
#include "Ysn_Fsn.h"
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
//RHS %21
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
	printf ("#Kernel=21\t#Threads=%d\t%.20e\n", omp_get_num_threads (),
		T / 1e9 / n);
#else
	printf ("#Kernel=21\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
      }
  }
#endif
#pragma omp barrier
//LC %13
#ifdef INSTRUMENT
  {
#pragma omp barrier
    time_snap_t time;
    if (me == 0)
      {
	time_snap_start (&time);
      }
#endif
    int bs_LC_lij, limit_LC_lij;
    for (int jj = first, bs_LC_lij =
	 imin (B_LC_lij, last - first + 1), limit_LC_lij =
	 imax (first, last + 1 - B_LC_lij); jj <= last;
	 bs_LC_lij = last + 1 - jj, limit_LC_lij = last)
      {
	for (; jj <= limit_LC_lij; jj += bs_LC_lij)
	  {
#pragma nounroll_and_jam
	    for (int l = 0; l < 4; ++l)
	      {
#pragma ivdep
		for (int j = jj; j < bs_LC_lij + jj; ++j)
		  {
		    Y[l][j] = A[l][0] * F[0][j];
		  }
#pragma nounroll_and_jam
		for (int i = 1; i < 4; ++i)
		  {
#pragma ivdep
		    for (int j = jj; j < bs_LC_lij + jj; ++j)
		      {
			Y[l][j] += A[l][i] * F[i][j];
		      }
		  }
#pragma ivdep
		for (int j = jj; j < bs_LC_lij + jj; ++j)
		  {
		    Y[l][j] = Y[l][j] * h + y[j];
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
#pragma omp barrier
  for (int k = 0; k < 5; ++k)
    {
//RHS %21
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
	    printf ("#Kernel=21\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=21\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
#pragma omp barrier
//LC %13
#ifdef INSTRUMENT
      {
#pragma omp barrier
	time_snap_t time;
	if (me == 0)
	  {
	    time_snap_start (&time);
	  }
#endif
	int bs_LC_lij, limit_LC_lij;
	for (int jj = first, bs_LC_lij =
	     imin (B_LC_lij, last - first + 1), limit_LC_lij =
	     imax (first, last + 1 - B_LC_lij); jj <= last;
	     bs_LC_lij = last + 1 - jj, limit_LC_lij = last)
	  {
	    for (; jj <= limit_LC_lij; jj += bs_LC_lij)
	      {
#pragma nounroll_and_jam
		for (int l = 0; l < 4; ++l)
		  {
#pragma ivdep
		    for (int j = jj; j < bs_LC_lij + jj; ++j)
		      {
			Y[l][j] = A[l][0] * F[0][j];
		      }
#pragma nounroll_and_jam
		    for (int i = 1; i < 4; ++i)
		      {
#pragma ivdep
			for (int j = jj; j < bs_LC_lij + jj; ++j)
			  {
			    Y[l][j] += A[l][i] * F[i][j];
			  }
		      }
#pragma ivdep
		    for (int j = jj; j < bs_LC_lij + jj; ++j)
		      {
			Y[l][j] = Y[l][j] * h + y[j];
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
	    printf ("#Kernel=13\t#Threads=%d\t%.20e\n",
		    omp_get_num_threads (), T / 1e9 / n);
#else
	    printf ("#Kernel=13\t#Threads=1\t%.20e\n", T / 1e9 / n);
#endif
	  }
      }
#endif
#pragma omp barrier
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
	 imin (B_RHS_Approx_Update_ij, last - first + 1),
	 limit_RHS_Approx_Update_ij =
	 imax (first, last + 1 - B_RHS_Approx_Update_ij); jj <= last;
	 bs_RHS_Approx_Update_ij = last + 1 - jj, limit_RHS_Approx_Update_ij =
	 last)
      {
	for (; jj <= limit_RHS_Approx_Update_ij;
	     jj += bs_RHS_Approx_Update_ij)
	  {
#pragma nounroll_and_jam
	    for (int i = 0; i < 4; ++i)
	      {
#pragma ivdep
		for (int j = jj; j < bs_RHS_Approx_Update_ij + jj; ++j)
		  {
		    y[j] +=
		      h * b[i] * (eval_component (j, t + c[i] * h, Y[i]));
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
