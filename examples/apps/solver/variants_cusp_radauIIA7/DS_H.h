#pragma once

double *y;
double **Y;
double *dy;
double **F;

static void
allocate_data_structures ()
{
  y = alloc1d (n);
  Y = alloc2d (s, n);
  dy = alloc1d (n);
  F = alloc2d (s, n);
}

static void
free_data_structures ()
{
  free1d (y);
  free2d (Y);
  free1d (dy);
  free2d (F);
}
