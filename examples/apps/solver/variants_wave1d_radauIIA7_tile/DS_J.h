#pragma once

double **F;
double *y;
double *Y;
double **Fprev;

static void
allocate_data_structures ()
{
  F = alloc2d (s, n);
  y = alloc1d (n);
  Y = alloc1d (n);
  Fprev = alloc2d (s, n);
}

static void
free_data_structures ()
{
  free2d (F);
  free1d (y);
  free1d (Y);
  free2d (Fprev);
}
