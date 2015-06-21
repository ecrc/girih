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
  printf("Global domain    size: %lu    nx:%d    ny:%d    nz:%d\n", l_stencils, sx, sy, sz); \
  printf("Number of time steps: %d\n", Nt); \
  printf("Number of tests: %d\n", TESTS); \
  printf("RANK0 MStencil/s  MAX: %f\n", perf); \
}


