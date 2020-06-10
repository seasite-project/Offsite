double **Yprev;
double **Ycur;
double **F;
double *y;
double *dy;
double *f;
double *Fs;
double *Fn;

static void allocate_data_structures()
{
  Yprev = alloc2d(s, n);
  Ycur = alloc2d(s, n);
  y = alloc1d(n);
  F = alloc2d(s, n);
  dy = alloc1d(n);
  f = alloc1d(1);
  Fs = alloc1d(s);
  Fn = alloc1d(n);
}

static void free_data_structures()
{
  free2d(Yprev);
  free2d(Ycur);
  free1d(y);
  free2d(F);
  free1d(dy);
  free1d(f);
  free1d(Fs);
  free1d(Fn);
}
