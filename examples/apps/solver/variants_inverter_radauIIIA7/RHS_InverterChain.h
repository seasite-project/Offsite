#pragma once

#define PROBLEM_NAME "InverterChain"
#define PROBLEM_ID 1

static inline void
eval_range (int first, int last, double t, const double *y, double *f)
{
#pragma ivdep
  for (int j = first; j <= last; j++)
    {
      double U_previous =
	((j ==
	  0) ? 1.00000000000000 * 5.00000000000000 * sin (2 * M_PI * t *
							  0.0125000000000000 -
							  M_PI / 2) * 0.5 +
	 (1.00000000000000 * 5.00000000000000 * 0.5) : y[j - 1]);
      double U_j = y[j];

      double max1 = max (U_previous - 1.00000000000000, 0.0);
      double max2 = max (U_previous - U_j - 1.00000000000000, 0.0);
      double ids = 0.000200000000000000 * (max1 * max1 - max2 * max2);

      f[j] =
	((5.00000000000000 - U_j) / 5000.00000000000 -
	 ids) / 0.00100000000000000;
    }
}

static inline double
eval_component (int j, double t, const double *y)
{
  double U_previous =
    ((j ==
      0) ? 1.00000000000000 * 5.00000000000000 * sin (2 * M_PI * t *
						      0.0125000000000000 -
						      M_PI / 2) * 0.5 +
     (1.00000000000000 * 5.00000000000000 * 0.5) : y[j - 1]);
  double U_j = y[j];

  double max1 = max (U_previous - 1.00000000000000, 0.0);
  double max2 = max (U_previous - U_j - 1.00000000000000, 0.0);
  double ids = 0.000200000000000000 * (max1 * max1 - max2 * max2);

  return ((5.00000000000000 - U_j) / 5000.00000000000 -
	  ids) / 0.00100000000000000;
}

static inline void
initial_values (const double t, double *y)
{
  for (int j = 0; j < n; j++)
    y[j] = ((j + 1) % 2) * 1.00000000000000 * 5.00000000000000;
}
