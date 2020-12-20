#pragma once

double **F;
double *y;
double **Y;

static void
allocate_data_structures ()
{
  F = alloc2d (s, n);
  y = alloc1d (n);
  Y = alloc2d (s, n);
}

static void
free_data_structures ()
{
  free2d (F);
  free1d (y);
  free2d (Y);
}
