#pragma once

#define PROBLEM_NAME "Heat2D"
#define PROBLEM_ID 4

static inline void
eval_range (int first, int last, double t, const double *y, double *f)
{
  int j = first;
  int ri = j / g;
  int rj = j / g;
  int ci = j % g;		// column of j
  int si = j - ci;		// first column in row of j
  int ei = si + g - 1;		// last column in row of j

  (void) t;

  if (ri == 0)
    {
      // first row is top row (maybe incomplete)
      if (ci == 0)
	{
	  f[j] = ((g - 1.0) * (g - 1.0)) * 2.0 * (-2.0 * y[0] + y[1] + y[g]);
	  j++;
	}

      int ti = imin (ei - 1, j);
#pragma ivdep
      for (; j <= ti; j++)
	f[j] =
	  ((g - 1.0) * (g - 1.0)) * (y[j - 1] - 4.0 * y[j] + y[j + 1] +
				     2.0 * y[j + g]);

      if (j > last)
	return;

      f[j] =
	((g - 1.0) * (g - 1.0)) * 2.0 * (y[g - 2] - 2.0 * y[g - 1] +
					 y[g + g - 1]);
      j++;

      ri++;
      ci = 0;

      if (ri > rj)
	return;

      si += g;
      ei += g;
    }

  if (ri < g - 1)
    {
      if (ci == 0)
	{
	  f[j] =
	    ((g - 1.0) * (g - 1.0)) * (y[j - g] + y[j + g] + 2.0 * y[j + 1] -
				       4.0 * y[j]);
	  j++;
	}

      int ti = imin (ei - 1, j);
#pragma ivdep
      for (; j <= ti; j++)
	f[j] =
	  ((g - 1.0) * (g - 1.0)) * (y[j - g] + y[j - 1] - 4.0 * y[j] +
				     y[j + 1] + y[j + g]);

      if (j > last)
	return;

      f[j] =
	((g - 1.0) * (g - 1.0)) * (y[j - g] + 2.0 * y[j - 1] - 4.0 * y[j] +
				   y[j + g]);
      j++;

      ri++;
      ci = 0;

      if (ri > rj)
	return;

      si += g;
      ei += g;
    }

  if (ri < g - 1)
    {
// complete inner rows
#pragma ivdep
      for (; ri < rj; ri++, si += g, ei += g)
	{
	  f[si] =
	    ((g - 1.0) * (g - 1.0)) * (y[j - g] + y[j + g] + 2.0 * y[j + 1] -
				       4.0 * y[j]);
	  j++;

	  for (; j < ei; j++)
	    f[j] =
	      ((g - 1.0) * (g - 1.0)) * (y[j - g] + y[j - 1] - 4.0 * y[j] +
					 y[j + 1] + y[j + g]);

	  f[ei] =
	    ((g - 1.0) * (g - 1.0)) * (y[j - g] + 2.0 * y[j - 1] -
				       4.0 * y[j] + y[j + g]);
	  j++;
	}
    }

  if (ri < g - 1)		// now, ri == rj
    {
      // last row is an inner row (maybe incomplete)
      f[si] =
	((g - 1.0) * (g - 1.0)) * (y[j - g] + y[j + g] + 2.0 * y[j + 1] -
				   4.0 * y[j]);
      j++;

      int ti = imin (ei - 1, j);
#pragma ivdep
      for (; j <= ti; j++)
	f[j] =
	  ((g - 1.0) * (g - 1.0)) * (y[j - g] + y[j - 1] - 4.0 * y[j] +
				     y[j + 1] + y[j + g]);

      if (j > last)
	return;

      f[ei] =
	((g - 1.0) * (g - 1.0)) * (y[j - g] + 2.0 * y[j - 1] - 4.0 * y[j] +
				   y[j + g]);
    }
  else
    {
      // last row is the bottom row (maybe incomplete)

      if (j == si)
	{
	  f[si] =
	    ((g - 1.0) * (g - 1.0)) * 2.0 * (y[n - 2 * g] - 2.0 * y[n - g] +
					     y[n - g + 1]);
	  j++;
	}

      int ti = imin (ei - 1, j);
#pragma ivdep
      for (; j <= ti; j++)
	f[j] =
	  ((g - 1.0) * (g - 1.0)) * (2.0 * y[j - g] + y[j - 1] - 4.0 * y[j] +
				     y[j + 1]);

      if (j > last)
	return;

      f[ei] =
	((g - 1.0) * (g - 1.0)) * 2.0 * (y[n - 1 - g] + y[n - 2] -
					 2.0 * y[n - 1]);
    }
}

static inline double
eval_component (int j, double t, const double *y)
{
  const int ri = j / g;
  const int ci = j % g;

  (void) t;

  if (ri > 0)
    {				// not the top row
      if (ri < g - 1)
	{			// not the bottom row
	  if (ci > 0)
	    {			// not the left column
	      if (ci < g - 1)
		return ((g - 1.0) * (g - 1.0)) * (y[j - g] + y[j - 1] -
						  4.0 * y[j] + y[j + 1] +
						  y[j + g]);
	      else		// ci == g - 1
		return ((g - 1.0) * (g - 1.0)) * (y[j - g] + 2.0 * y[j - 1] -
						  4.0 * y[j] + y[j + g]);
	    }
	  else			// ci == 0
	    return ((g - 1.0) * (g - 1.0)) * (y[j - g] + y[j + g] +
					      2.0 * y[j + 1] - 4.0 * y[j]);
	}
      else			// ri == g - 1
	{			// bottom row
	  if (ci > 0)
	    {			// not the bottom left corner
	      if (ci < g - 1)
		return ((g - 1.0) * (g - 1.0)) * (2.0 * y[j - g] + y[j - 1] -
						  4.0 * y[j] + y[j + 1]);
	      else		// ci == g - 1
		return ((g - 1.0) * (g - 1.0)) * 2.0 * (y[n - 1 - g] +
							y[n - 2] - 2.0 * y[n -
									   1]);
	    }
	  else			// ci == 0
	    return ((g - 1.0) * (g - 1.0)) * 2.0 * (y[n - 2 * g] -
						    2.0 * y[n - g] + y[n - g +
								       1]);
	}
    }
  else				// ri == 0
    {				// top row
      if (ci > 0)
	{			// not the top left corner
	  if (ci < g - 1)
	    return ((g - 1.0) * (g - 1.0)) * (y[j - 1] - 4.0 * y[j] +
					      y[j + 1] + 2.0 * y[j + g]);
	  else			// ci == g - 1
	    return ((g - 1.0) * (g - 1.0)) * 2.0 * (y[g - 2] -
						    2.0 * y[g - 1] + y[g + g -
								       1]);
	}
      else			// ci == 0
	return ((g - 1.0) * (g - 1.0)) * 2.0 * (-2.0 * y[0] + y[1] + y[g]);
    }
}

static inline void
initial_values (const double t, double *y)
{
  for (int i = 0; i < g; i++)
    for (int j = 0; j < g; j++)
      y[i * g + j] = 0.5 + (double) j / (double) (g - 1);
}
