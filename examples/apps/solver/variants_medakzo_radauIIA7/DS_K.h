#pragma once

double *y;
double *dy;
double *Y;
double **Fprev;
double **F;

static void
allocate_data_structures ()
{
  y = alloc1d (n);
  dy = alloc1d (n);
  Y = alloc1d (n);
  Fprev = alloc2d (s, n);
  F = alloc2d (s, n);
}

static void
free_data_structures ()
{
  free1d (y);
  free1d (dy);
  free1d (Y);
  free2d (Fprev);
  free2d (F);
}
