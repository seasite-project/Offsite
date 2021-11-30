#pragma once

double **F;
double *y;
double **Y;
double *Yn;
double tmp;

static void
allocate_data_structures ()
{
  F = alloc2d (s, n);
  y = alloc1d (n);
  Y = alloc2d (s, n);
  Yn = alloc1d (s);
}

static void
free_data_structures ()
{
  free2d (F);
  free1d (y);
  free2d (Y);
  free1d (Yn);
}
