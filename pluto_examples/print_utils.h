#define PRINT_RESULTS(r, coef) { \
  int sx = Nx-2*r; \
  int sy = Ny-2*r; \
  int sz = Nz-2*r; \
  unsigned long int l_stencils = 1UL*sx*sy*sz; \
  double perf = (l_stencils*Nt*1e-6)/(min_tdiff); \
  printf("Time stepper name: PLUTO\n"); \
  printf("Stencil Kernel name: star\n"); \
  printf("Stencil Kernel semi-bandwidth: %d\n", (r)); \
  printf("Stencil Kernel coefficients: %s\n", (coef)); \
  printf("Global domain    size:%lu    nx:%d    ny:%d    nz:%d\n", l_stencils, sx, sy, sz); \
  printf("Number of time steps: %d\n", Nt); \
  printf("Number of tests: %d\n", TESTS); \
  printf("OpenMP Threads: %d\n", num_threads); \
  int ti=0; \
  while(tile_size[ti]!=-1) {\
     printf("PLUTO tile size of loop %d: %d\n", ti+1, tile_size[ti]); \
    ti++; \
  } \
  printf("RANK0 MStencil/s  MAX: %f\n", perf); \
}

//  char *t_str=NULL; \
  size_t t_len=0; \
  int idx=0; \
  FILE * fp = fopen("tile.sizes", "r"); \
  if(fp!=NULL){ \
      while ( getline(&t_str, &t_len, fp) != -1) { \
          printf("PLUTO tile size of loop %d: %s", idx++, t_str); \
      } \
      if (t_str) free(t_str); \
      fclose(fp); \
  } 
