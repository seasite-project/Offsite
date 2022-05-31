#pragma once

double *y;
double *dy;
double **F;
double **Y;
double **Yprev;
double **Ycur;
double *tmp;
double f;
double *Fs;
double *Fn;
double tmp_s;

static void
allocate_data_structures ()
{
  y = alloc1d (n);
  dy = alloc1d (n);
  F = alloc2d (s, n);
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
  free1d (y);
  free1d (dy);
  free2d (F);
  free2d (Y);
  free2d (Yprev);
  free2d (Ycur);
  free1d (tmp);
  free1d (Fs);
  free1d (Fn);
}
