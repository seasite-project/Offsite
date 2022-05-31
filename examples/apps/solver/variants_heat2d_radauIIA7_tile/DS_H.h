#pragma once

double *y;
double **Y;
double dy;
double **F;
double *Yn;
double tmp;

static void
allocate_data_structures ()
{
  y = alloc1d (n);
  Y = alloc2d (s, n);
  F = alloc2d (s, n);
  Yn = alloc1d (s);
}

static void
free_data_structures ()
{
  free1d (y);
  free2d (Y);
  free2d (F);
  free1d (Yn);
}
