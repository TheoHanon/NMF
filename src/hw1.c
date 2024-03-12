#include "../header/RK4.h"

#include "../header/fd_scheme.h"

#include "../header/struct.h"


int main(int argc, char * argv[]) {
  (void) argc;
  const int N = atoi(argv[1]);
  const SchemeType my_scheme_type = get_scheme_type(argv[2]);
  const GridType my_grid_type = get_grid_type(argv[3]);
  const InitialType my_initial_type = get_initial_type(argv[4]);

  char * file_name1 = argv[5];
  char * file_name2 = argv[6];
  char * file_name3 = argv[7];
  char * file_name4 = argv[8];

  const double L = 1.0;
  const double dx = L / N;
  const double sigma = L / 16.0;
  const double a = 3.0 / 5.0;
  const double c = 1;
  const double beta = 0.5;
  const double dt = beta * dx / c;
  const double mu = 0.0;
  const int N_step = 5000;
  const double U = 1.0;
  const int p = 12;

  Gaussian * g = init_gaussian(U, mu, sigma, L);
  Scheme * s = init_scheme(dt, dx, N, my_scheme_type);
  Grid * my_grid = init_grid(dx, N, my_grid_type, a, L);
  Parameter * my_param = init_parameter(N, L, dx, a, dt, beta, c, N_step, U, my_grid_type, my_scheme_type);

  printf("==== Simulation parameters ====\n");
  printf("N = %d\n", N);
  printf("L = %.3f\n", L);
  printf("sigma = %.3f\n", g -> sigma);
  printf("dx = %.3f\n", dx);
  printf("dt = %.3f\n", dt);
  printf("c = %.3f\n", c);
  printf("N_step = %d\n", N_step);
  printf("U = %.3f\n", U);
  printf("mu = %.3f\n", g -> mu);
  printf("sigma = %.3f\n", g -> sigma);
  printf("CFL = %.3f\n", beta);

  if (my_grid -> type == NONUNIF) {
    printf("Grid type = NON_UNIF\n");
  } else {
    printf("Grid type = UNIF\n");
  }

  if (s -> type == ED) {
    printf("Scheme type = ED\n");
  } else if (s -> type == E2) {
    printf("Scheme type = E2\n");
  } else if (s -> type == E4) {
    printf("Scheme type = E4\n");
  } else {
    printf("Scheme type = I4\n");
  }

  if (my_initial_type == GAUSSIAN) {
    printf("Initial condition = GAUSSIAN\n");
  } else if (my_initial_type == WAVEPACKET) {
    printf("Initial condition = WAVEPACKET\n");
  } else {
    printf("Invalid initial condition\n");
    exit(1);
  }
  printf("================================\n");

  double * u0 = (double * ) malloc(N * sizeof(double));
  // Dynamically allocate memory for the array
  // Initialize the array
  save_grid(my_grid, file_name3);
  save_parameter(my_param, file_name4);

  if (my_initial_type == GAUSSIAN) {
    initGaussianArray( & u0, N, my_grid, g);
  } else if (my_initial_type == WAVEPACKET) {
    initWavePacketArray( & u0, N, p, g, my_grid);
  } else {
    printf("Invalid initial condition\n");
    exit(1);
  }
  convection_solve(file_name1, file_name2, & u0, N_step, c, my_grid, s, g);

  free(u0);
  free_grid(my_grid);
  free_scheme(s);
  free_parameter(my_param);
  free(g);

  return 0;
}

