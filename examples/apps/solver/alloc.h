#define ALIGNMENT 32

#define ALIGN(X) ((unsigned long) ((((X) + ALIGNMENT - 1)/ALIGNMENT) * ALIGNMENT) % 256 < ALIGNMENT ? ((((X) + ALIGNMENT - 1)/ALIGNMENT) * ALIGNMENT) + ALIGNMENT : ((((X) + ALIGNMENT - 1)/ALIGNMENT) * ALIGNMENT))

static double *alloc1d(size_t a)
{
  return (double *) aligned_alloc(ALIGNMENT, a * sizeof(double));
}

static double **alloc2d(size_t a, size_t b)
{
  size_t i, row_size, row_count;
  double **x ;

  row_size = ALIGN(b * sizeof(double));
  row_count = row_size / sizeof(double);
  
  x = (double **) aligned_alloc(ALIGNMENT, a * sizeof(double *));
  x[0] = (double *) aligned_alloc(ALIGNMENT, a * row_size);
  for (i = 1; i < a; i++)
    x[i] = x[0] + i * row_count;

  return x;
}

static double ***alloc3d(size_t a, size_t b, size_t c)
{
  size_t i, row_size, row_count;
  double ***x ;

  assert((ALIGNMENT % sizeof(double)) == 0);

  row_size = ALIGN(c * sizeof(double));
  row_count = row_size / sizeof(double);
  
  x = (double ***) aligned_alloc(ALIGNMENT, a * sizeof(double **));
  x[0] = (double **) aligned_alloc(ALIGNMENT, a * b * sizeof(double *));
  x[0][0] = (double *) aligned_alloc(ALIGNMENT, a * b * row_size);

  for (i = 1; i < a; i++)
    x[i] = x[0] + i * b;

  for (i = 1; i < a * b; i++)
      x[0][i] = x[0][0] + i * row_count;

  return x;
}

static void free1d(double *p)
{
  free((void *)p);
}

static void free2d(double **p)
{
  free((void *)p[0]);
  free((void *)p);
}

static void free3d(double ***p)
{
  free((void *)p[0][0]);
  free((void *)p[0]);
  free((void *)p);
}

