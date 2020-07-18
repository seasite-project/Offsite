#pragma once

double *y;
double *dy;
double **Yprev;
double **Ycur;
double *f;
double *Fs;
double *Fn;

static void
allocate_data_structures ()
{
  y = alloc1d (n);
  dy = alloc1d (n);
  Yprev = alloc2d (s, n);
  Ycur = alloc2d (s, n);
  f = alloc1d (1);
  Fs = alloc1d (s);
  Fn = alloc1d (n);
}

static void
free_data_structures ()
{
  free1d (y);
  free1d (dy);
  free2d (Yprev);
  free2d (Ycur);
  free1d (f);
  free1d (Fs);
  free1d (Fn);
}
