#pragma once

#define PROBLEM_NAME "Medakzo"
#define PROBLEM_ID 3

static inline void
eval_range (int first, int last, double t, const double *y, double *f)
{
  int j = first;

  if (first < g && last >= g)
    {				// Both PDEs
      if (first == 0)
	{
	  double phi = 0.0;
	  if (t <= 5.0)
	    phi = 2.0;

	  double tmp = (1.0 / g - 1.0) * (1.0 / g - 1.0) / 4.00000000000000;
	  double alpha = 2.0 * (1.0 / g - 1.0) * tmp / 4.00000000000000;
	  double beta = tmp * tmp;

	  f[j] =
	    alpha * (y[1] - phi) / (2.0 * 1.0 / g) + (phi - 2.0 * y[0] +
						      y[1]) * beta / (1.0 /
								      g *
								      1.0 /
								      g) -
	    100.000000000000 * y[0] * y[g];
	  ++j;
	}

#pragma ivdep
      for (; j < g - 1; ++j)
	{
	  double zeta = ((j + 1) / 2) * 1.0 / g;
	  double tmp = (zeta - 1.0) * (zeta - 1.0) / 4.00000000000000;
	  double alpha = 2.0 * (zeta - 1.0) * tmp / 4.00000000000000;
	  double beta = tmp * tmp;
	  f[j] =
	    (alpha * (y[j + 1] - y[j - 1]) / (2.0 * 1.0 / g) +
	     (y[j - 1] - 2.0 * y[j] +
	      y[j + 1]) * beta / (1.0 / g * 1.0 / g)) -
	    100.000000000000 * y[j] * y[j + g];
	}

      f[j] = -100.000000000000 * y[g - 1] * y[n - 1];
      ++j;

      f[j] = -100.000000000000 * y[0] * y[g];
      ++j;

      if (last < n - 1)
	{
#pragma ivdep
	  for (; j <= last; ++j)
	    f[j] = -100.000000000000 * y[j] * y[j - g];
	}
      else if (last == n - 1)
	{
#pragma ivdep
	  for (; j < n - 1; ++j)
	    f[j] = -100.000000000000 * y[j] * y[j - g];
	  f[j] = (-100.000000000000 * y[g - 1] * y[n - 1]);
	}
    }
  if (last < g)			// Only 1.PDE
    {
      if (first == 0)
	{
	  double phi = 0.0;
	  if (t <= 5.0)
	    phi = 2.0;

	  double tmp = (1.0 / g - 1.0) * (1.0 / g - 1.0) / 4.00000000000000;
	  double alpha = 2.0 * (1.0 / g - 1.0) * tmp / 4.00000000000000;
	  double beta = tmp * tmp;

	  f[j] =
	    alpha * (y[1] - phi) / (2.0 * 1.0 / g) + (phi - 2.0 * y[0] +
						      y[1]) * beta / (1.0 /
								      g *
								      1.0 /
								      g) -
	    100.000000000000 * y[0] * y[g];
	  ++j;
	}

      if (last < g - 1)
	{
#pragma ivdep
	  for (; j <= last; ++j)
	    {
	      double zeta = ((j + 1) / 2) * 1.0 / g;
	      double tmp = (zeta - 1.0) * (zeta - 1.0) / 4.00000000000000;
	      double alpha = 2.0 * (zeta - 1.0) * tmp / 4.00000000000000;
	      double beta = tmp * tmp;
	      f[j] =
		(alpha * (y[j + 1] - y[j - 1]) / (2.0 * 1.0 / g) +
		 (y[j - 1] - 2.0 * y[j] +
		  y[j + 1]) * beta / (1.0 / g * 1.0 / g)) -
		100.000000000000 * y[j] * y[j + g];
	    }
	}
      else if (last == g - 1)
	{
#pragma ivdep
	  for (; j < g - 1; ++j)
	    {
	      double zeta = ((j + 1) / 2) * 1.0 / g;
	      double tmp = (zeta - 1.0) * (zeta - 1.0) / 4.00000000000000;
	      double alpha = 2.0 * (zeta - 1.0) * tmp / 4.00000000000000;
	      double beta = tmp * tmp;
	      f[j] =
		(alpha * (y[j + 1] - y[j - 1]) / (2.0 * 1.0 / g) +
		 (y[j - 1] - 2.0 * y[j] +
		  y[j + 1]) * beta / (1.0 / g * 1.0 / g)) -
		100.000000000000 * y[j] * y[j + g];
	    }

	  f[j] = -100.000000000000 * y[g - 1] * y[n - 1];
	}
    }
  else if (first >= g)		// 2.PDE
    {
      if (first == g)
	{
	  f[j] = -100.000000000000 * y[0] * y[g];
	  ++j;
	}

      if (last < n - 1)
	{
#pragma ivdep
	  for (; j <= last; ++j)
	    f[j] = -100.000000000000 * y[j] * y[j - g];
	}
      else
	{
#pragma ivdep
	  for (; j < n - 1; ++j)
	    f[j] = -100.000000000000 * y[j] * y[j - g];

	  f[j] = (-100.000000000000 * y[g - 1] * y[n - 1]);
	}
    }
}

static inline double
eval_component (int j, double t, const double *y)
{
  if (j > 0 && j < g - 1)
    {
      double zeta = ((j + 1) / 2) * 1.0 / g;
      double tmp = (zeta - 1.0) * (zeta - 1.0) / 4.00000000000000;
      double alpha = 2.0 * (zeta - 1.0) * tmp / 4.00000000000000;
      double beta = tmp * tmp;

      return (alpha * (y[j + 1] - y[j - 1]) / (2.0 * 1.0 / g) +
	      (y[j - 1] - 2.0 * y[j] +
	       y[j + 1]) * beta / (1.0 / g * 1.0 / g)) -
	100.000000000000 * y[j] * y[j + g];
    }
  else if (j == 0)
    {
      double phi = 0.0;
      if (t <= 5.0)
	phi = 2.0;

      double tmp = (1.0 / g - 1.0) * (1.0 / g - 1.0) / 4.00000000000000;
      double alpha = 2.0 * (1.0 / g - 1.0) * tmp / 4.00000000000000;
      double beta = tmp * tmp;

      return alpha * (y[1] - phi) / (2.0 * 1.0 / g) + (phi - 2.0 * y[0] +
						       y[1]) * beta / (1.0 /
								       g *
								       1.0 /
								       g) -
	100.000000000000 * y[0] * y[g];
    }
  else if (j == g - 1)
    return -100.000000000000 * y[g - 1] * y[n - 1];
  else if (j == g)
    return -100.000000000000 * y[0] * y[g];
  else if (j < 2 * g - 1)
    return -100.000000000000 * y[j] * y[j - g];
  else if (j == 2 * g - 1)
    return (-100.000000000000 * y[g - 1] * y[n - 1]);
  return -1;
}

static inline void
initial_values (const double t, double *y)
{
  for (int i = 1; i <= n / 2; i++)
    {
      y[i - 1] = 0.0;
      y[g + i - 1] = 1.00000000000000;
    }
}
