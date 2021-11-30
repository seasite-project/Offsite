#pragma once

double **F;
double *y;
double **Y;
double **Yprev;
double **Ycur;
double *tmp;
double f;
double *Fs;
double *Fn;

static void
allocate_data_structures ()
{
  F = alloc2d (s, n);
  y = alloc1d (n);
  Y = alloc2d (s, n);
  Yprev = alloc2d (s, n);
  Ycur = alloc2d (s, n);
  tmp = alloc1d (s);
  Fs = alloc1d (s);
  Fn = alloc1d (n);
}

static void
free_data_structures ()
{
  free2d (F);
  free1d (y);
  free2d (Y);
  free2d (Yprev);
  free2d (Ycur);
  free1d (tmp);
  free1d (Fs);
  free1d (Fn);
}
