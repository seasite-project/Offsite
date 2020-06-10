#ifndef TIMESNAP_H
#define TIMESNAP_H

///////////////////////////////////////////////////////////////////////

#if defined(TIMESNAP_PAPI)

#include <stdint.h>
#include <papi.h>

typedef long long time_snap_t;

static inline void init_time_snap()
{
  int events[] = {PAPI_TOT_CYC};
  PAPI_start_counters(events, 1);
}

static inline long long papi_timesnap()
{
  long long values[1];
  PAPI_read_counters(values, 1);
  return values[0];
}

static inline void time_snap_start(time_snap_t *ts)
{
  *ts = papi_timesnap();
}

static inline uint64_t time_snap_stop(const time_snap_t *ts1)
{
  time_snap_t ts2 = papi_timesnap();
  return ((ts2 < 0) ? 0 : (uint64_t) ts2);
}

#elif defined(TIMESNAP_CLOCK_GETTIME)

#include <time.h>

#ifdef _POSIX_MONOTONIC_CLOCK

typedef struct timespec time_snap_t;

static inline void init_time_snap()
{
  time_snap_t ts;
  int err;

  err = clock_getres(CLOCK_MONOTONIC, &ts);

  printf("Timer resolution is: %ld s, %ld ns (return code %d)\n", ts.tv_sec,
         ts.tv_nsec, err);

  if (err != 0)
    libsolve_abort(1);
}

static inline void time_snap_start(time_snap_t *ts)
{
  clock_gettime(CLOCK_MONOTONIC, ts);
}

static inline uint64_t time_snap_stop(const time_snap_t *ts1)
{
  time_snap_t ts2;
  clock_gettime(CLOCK_MONOTONIC, &ts2);

  return (uint64_t)(ts2.tv_sec - ts1->tv_sec) * 1000000000 +
         (uint64_t) ts2.tv_nsec - (uint64_t) ts1->tv_nsec;
}

#else

#error clock_gettime(CLOCK_MONOTONIC, ...) is not available.

#endif

#else /* TIMESNAP_GETTIMEOFDAY */

#include <sys/time.h>
#include <stdint.h>

typedef struct timeval time_snap_t;

#define init_time_snap()

static inline void time_snap_start(time_snap_t *ts)
{
  gettimeofday(ts, NULL);
}

static inline uint64_t time_snap_stop(const time_snap_t *ts1)
{
  time_snap_t ts2;
  gettimeofday(&ts2, NULL);

  return ((uint64_t)(ts2.tv_sec - ts1->tv_sec) * 1000000 +
          (uint64_t) ts2.tv_usec - (uint64_t) ts1->tv_usec) *
         1000;
}

#endif

///////////////////////////////////////////////////////////////////////

#endif
