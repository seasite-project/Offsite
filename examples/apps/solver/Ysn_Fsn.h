double **Y;
double **F;
double *y;
double *dy;

static void allocate_data_structures()
{
  Y = alloc2d(s, n);
  y = alloc1d(n);
  F = alloc2d(s, n);
  dy = alloc1d(n);
}

static void free_data_structures()
{
  free2d(Y);
  free1d(y);
  free2d(F);
  free1d(dy);
}
