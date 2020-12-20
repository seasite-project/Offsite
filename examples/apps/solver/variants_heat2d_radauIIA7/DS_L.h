#pragma once

double *Y;
double **Fprev;
double *y;
double **F;

static void
allocate_data_structures ()
{
  Y = alloc1d (n);
  Fprev = alloc2d (s, n);
  y = alloc1d (n);
  F = alloc2d (s, n);
}

static void
free_data_structures ()
{
  free1d (Y);
  free2d (Fprev);
  free1d (y);
  free2d (F);
}
