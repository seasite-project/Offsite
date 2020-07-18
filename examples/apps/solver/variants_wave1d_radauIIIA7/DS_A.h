#pragma once

double *y;
double *dy;
double **F;
double **Y;

static void
allocate_data_structures ()
{
  y = alloc1d (n);
  dy = alloc1d (n);
  F = alloc2d (s, n);
  Y = alloc2d (s, n);
}

static void
free_data_structures ()
{
  free1d (y);
  free1d (dy);
  free2d (F);
  free2d (Y);
}
