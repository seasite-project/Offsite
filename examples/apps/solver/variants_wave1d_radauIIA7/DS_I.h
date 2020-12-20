#pragma once

double *y;
double *dy;
double **F;
double *Y;
double **Fprev;

static void
allocate_data_structures ()
{
  y = alloc1d (n);
  dy = alloc1d (n);
  F = alloc2d (s, n);
  Y = alloc1d (n);
  Fprev = alloc2d (s, n);
}

static void
free_data_structures ()
{
  free1d (y);
  free1d (dy);
  free2d (F);
  free1d (Y);
  free2d (Fprev);
}
