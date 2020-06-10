double *F;
double **k;
double *y;
double *dy;

static void allocate_data_structures()
{
  F = alloc1d(n);
  k = alloc2d(s, n);
  y = alloc1d(n);
  dy = alloc1d(n);
}

static void free_data_structures()
{
  free1d(F);
  free2d(k);
  free1d(y);
  free1d(dy);
}
