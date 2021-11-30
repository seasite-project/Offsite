#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#ifdef _OPENMP
#include <omp.h>
#else
#include <sched.h>
#endif

#ifndef THREADS
#define THREADS 1
#endif

double H;

static inline double dmin(double a, double b)
{
  if (a < b)
    return a;
  else
    return b;
}

static inline double dmax(double a, double b)
{
  if (a < b)
    return b;
  else
    return a;
}

static inline int imin(int a, int b)
{
  if (a < b)
    return a;
  else
    return b;
}

static inline int imax(int a, int b)
{
  if (a < b)
    return b;
  else
    return a;
}

#include "alloc.h"

#ifndef VARIANT_HEADER
#error "Variant not defined!"
#else
#include VARIANT_HEADER
#endif

#include "progress.h"
#include "timesnap.h"

int main()
{
  int steps = 0, nthreads;

  double hmin = h0;
  double hmax = h0;
  double hsum = 0.0;

  time_snap_t ts;
  double T, tps, tpcs;

  H = te - t0;

#ifdef _OPENMP
#pragma omp parallel num_threads(THREADS)
  {
#pragma omp single
    nthreads = omp_get_num_threads();
    if (nthreads != THREADS)
    {
      fprintf(stderr, "Conflicting number of threads!\n");
      exit(1);
    }
  }
#else
  { // When OpenMP is not used, take care of CPU affinity by ourselves.
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(0, &mask);
    sched_setaffinity(0, sizeof(cpu_set_t), &mask);
  }
#endif

#if defined(BENCHMARK) || defined(INSTRUMENT)
  {
    if (n % THREADS != 0)
    {
      printf("@\t@\t@\t@\t@\t@\t@\n");
      fprintf(stderr, "Number of equations (%d) is not an integer multiple of the number of threads (%d).\n", n, THREADS);
      exit(1);
    }
  }
#endif

#ifndef INSTRUMENT
#ifdef BENCHMARK
  printf("%s\t%s\t%d\t%d\t%s\t%d\t", VARIANT_NAME, METHOD_NAME, s, m, PROBLEM_NAME, n);
#ifdef _OPENMP
  printf("openmp\t");
#else
  printf("seq\t");
#endif
  printf("%d\t", THREADS);
#else
  printf("Variant: %s\n", VARIANT_NAME);
  printf("Method:  %s\n", METHOD_NAME);
  printf("Problem: %s\n", PROBLEM_NAME);
  printf("\n");
  printf("n = %d\n", n);
  printf("p = %d\n", THREADS);
  printf("\n");
#endif
#else 
  printf("#Variant=%d\t#Method=%d\t#IVP=%d\t#n=%d\n", VARIANT_ID, METHOD_ID, PROBLEM_ID, n);
#endif

  if (THREADS > n)
  {
    printf("@\t@\t@\t@\t@\t@\t@\n");
    fprintf(stderr, "Error: More threads than equations!\n");
    exit(0);
  }

  allocate_data_structures();

  initial_values(t0, y);

  print_pitch_line();

  init_time_snap();

  time_snap_start(&ts);

#pragma omp parallel num_threads(THREADS)
  {
    int me, first, last;
    double t = t0;
    double h = h0;

#ifdef _OPENMP
    me = omp_get_thread_num();
#else
    me = 0;
#endif

    first = me * n / THREADS;
    last = (me + 1) * n / THREADS - 1;

    while (t < te)
    {

#pragma omp master
      {
        steps++;
      
        hmin = dmin(hmin, h);
        hmax = dmax(hmax, h);
        hsum += h;
      }

#ifndef BENCHMARK
#pragma omp barrier

      if (h < 1e-12)
      {
#pragma omp master
        printf("%d: h=%e too small!\n", me, h);
        fflush(stdout);
        break;
      }

#endif

      h = dmin(h, te - t);
        
#pragma omp master
      show_progress(t, h);

      timestep(me, first, last, t, h);

      t += h;
    }
  }

  T = time_snap_stop(&ts) / 1e9;
  tps = T / (double) steps;
  tpcs = tps / (double) n;

  show_completion();

#ifndef INSTRUMENT
#ifdef BENCHMARK
  printf("%d\t%.20e\t%.20e\t%.20e\t%e\t%e\t%e\n",
         steps, hmin, hmax, hsum / (double) steps, T, tps, tpcs);
#else
  printf("\n");
  printf("steps = %d\n", steps);
  printf("\n");
  printf("hmin  = %.20e\n", hmin);
  printf("hmax  = %.20e\n", hmax);
  printf("havg  = %.20e\n", hsum / (double) steps);
  printf("\n");
  printf("time  = %.3e s\n", T);
  printf("tps   = %.3e s\n", tps);
  printf("tpcs  = %.3e s\n", tpcs);
#endif
#endif

  free_data_structures();

  return 0;
}
