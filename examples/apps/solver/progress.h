#ifdef BENCHMARK

#define print_pitch_line()
#define show_progress(t, h)
#define show_completion()

#else

int percent = 0;

void print_pitch_line()
{
  int i;

  for (i = 1; i <= 100; i++)
  {
    if (i % 10 == 0)
      printf(",");
    else
      printf("_");
  }
  printf("\n");
}

void show_progress(double t, double h)
{
  int i;
  int percent_new = (int) ((t - t0) / H * 100.0);

  if (percent_new != percent)
  {
    double step_ratio = h/H;

    for (i = percent + 1; i <= percent_new; i++)
    {
      if (step_ratio >= 1e-4)
      {
        printf("`"); 
      }
      else if (step_ratio >= 1e-5)
        printf("+"); 
      else if (step_ratio >= 1e-6)
        printf("*"); 
      else if (step_ratio >= 1e-7)
        printf("#"); 
      else
        printf("!");
    }
      
    fflush(stdout);
    percent = percent_new;
  }
}

void show_completion()
{
  int i;

  for (i = percent + 1; i <= 100; i++)
    printf("`"); 
  printf("\n");
}

#endif
