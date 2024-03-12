#include "../header/struct.h"

// Parser functions

GridType get_grid_type(const char * type) {
  if (strcmp(type, "UNIF") == 0) {
    return UNIF;
  } else {
    return NONUNIF;
  }
}

SchemeType get_scheme_type(const char * type) {
  if (strcmp(type, "ED") == 0) {
    return ED;
  } else if (strcmp(type, "E2") == 0) {
    return E2;
  } else if (strcmp(type, "E4") == 0) {
    return E4;
  } else {
    return I4;
  }
}

InitialType get_initial_type(const char * type) {
  if (strcmp(type, "GAUSSIAN") == 0) {
    return GAUSSIAN;
  } else {
    return WAVEPACKET;
  }
}

// Structure related functions

Gaussian * init_gaussian(double U, double mu, double sigma, double L) {
  Gaussian * g = malloc(sizeof(Gaussian));
  g -> U = U;
  g -> mu = mu;
  g -> sigma = sigma;
  g -> xmin = -L / 2.0;
  g -> xmax = L / 2.0;
  g -> L = L;
  return g;
}

Grid * init_grid(const double dx,
  const int N,
    const GridType type,
      const double a,
        const double L) {
  Grid * g = malloc(sizeof(Grid));
  g -> dx = dx;
  g -> N = N;
  g -> type = type;
  g -> a = a;
  g -> L = L;
  return g;
}

void free_grid(Grid * g) {
  free(g);
  return;
}

Parameter * init_parameter(int N, double L, double dx, double a, double dt, double beta, double c, int N_step, double U, GridType type, SchemeType scheme_type) {
  Parameter * p = malloc(sizeof(Parameter));
  p -> N = N;
  p -> L = L;
  p -> dx = dx;
  p -> a = a;
  p -> dt = dt;
  p -> beta = beta;
  p -> c = c;
  p -> N_step = N_step;
  p -> U = U;
  p -> type = type;
  p -> scheme_type = scheme_type;
  return p;
}

void free_parameter(Parameter * param) {
  free(param);
  return;
}

struct Scheme * init_scheme(const double dt,
  const double dx,
    const int N, SchemeType scheme_type) {

  struct Scheme * s = malloc(sizeof(Scheme));

  s -> dt = dt;
  s -> dx = dx;
  s -> N = N;
  s -> type = scheme_type;

  return s;
}

void free_scheme(struct Scheme * s) {
  free(s);
  return;
}

void save_grid(Grid * grid, char * file_name) {
  FILE * fptr = fopen(file_name, "w");
  int N = grid -> N;
  double L = grid -> L;
  double dx = grid -> dx;
  double x;
  for (int i = 0; i < N; i++) {
    x = -L / 2.0 + i * dx;
    fprintf(fptr, "%.5f\n", x);
  }
  fclose(fptr);
  return;
}

void save_parameter(Parameter * param, char * file_name) {
  FILE * fptr = fopen(file_name, "w");
  fprintf(fptr, "N=%d\n", param -> N);
  fprintf(fptr, "L=%.5f\n", param -> L);
  fprintf(fptr, "dx=%.5f\n", param -> dx);
  fprintf(fptr, "a=%.5f\n", param -> a);
  fprintf(fptr, "dt=%.5f\n", param -> dt);
  fprintf(fptr, "beta=%.5f\n", param -> beta);
  fprintf(fptr, "c=%.5f\n", param -> c);
  fprintf(fptr, "N_step=%d\n", param -> N_step);
  fprintf(fptr, "U=%.5f\n", param -> U);
  if (param -> type == NONUNIF) {
    fprintf(fptr, "Grid type=NON_UNIF\n");
  } else {
    fprintf(fptr, "Grid type=UNIF\n");
  }
  if (param -> scheme_type == ED) {
    fprintf(fptr, "Scheme type=ED\n");
  } else if (param -> scheme_type == E2) {
    fprintf(fptr, "Scheme type=E2\n");
  } else if (param -> scheme_type == E4) {
    fprintf(fptr, "Scheme type=E4\n");
  } else {
    fprintf(fptr, "Scheme type=I4\n");
  }
  fclose(fptr);
  return;
}

void save(FILE * fptr,
  const double * u,
    const double t,
      const int N) {
  fprintf(fptr, "%.12f", t);
  for (int j = 0; j < N; ++j) {
    fprintf(fptr, " %.12f", u[j]);
  }
  fprintf(fptr, "\n");
  return;
}

double gaussian(double x, struct Gaussian * g) {
  double sigma = g -> sigma;
  double mu = g -> mu;
  double U = g -> U;

  return U * exp(-pow(x - mu, 2) / (sigma * sigma));
}

void initGaussianArray(double ** arr, int n, Grid * grid, struct Gaussian * g) {
  double dx = grid -> dx;
  double L = grid -> L;
  double a = grid -> a;

  double x;
  for (int i = 0; i < n; i++) {
    x = -L / 2 + i * dx;
    if (grid -> type == NONUNIF) {

      ( * arr)[i] = gaussian(g_f(x, grid), g) * (1 - a * cos(2 * M_PI * x / L));
    } else {
      ( * arr)[i] = gaussian(x, g);
    }
  }
}

void initWavePacketArray(double ** arr, int n, int p, Gaussian * g, Grid * grid) {
  double x, k;
  double dx = grid -> dx;
  double L = grid -> L;
  double U = g -> U;
  double a = grid -> a;

  for (int i = 0; i < n; i++) {
    x = -L / 2 + i * dx;
    k = 2 * M_PI * p / L;
    if (grid -> type == NONUNIF) {
      ( * arr)[i] = U * cos(k * g_f(x, grid)) * gaussian(g_f(x, grid), g) * (1 - a * cos(2 * M_PI * x / L));
    } else {
      ( * arr)[i] = U * cos(k * x) * gaussian(x, g);
    }
  }
}

double gaussian_solution(const double x,
  const double t,
    const double c, struct Gaussian * g) {
  double L = g -> L;
  double xtild = fmod((x - c * t) + L / 2.0, L);
  xtild = (xtild < 0) ? xtild + L : xtild;
  xtild = xtild - L / 2.0;
  double sol = gaussian(xtild, g);
  return sol;
}

void diagnose(FILE * file, double ** u, double t,
  const int N,
    const double dx, struct Gaussian * g,
      const double c) {

  double I = 0.0;
  double E = 0.0;
  double R = 0.0;

  double sigma = g -> sigma;
  double U = g -> U;
  double x_min = g -> xmin;

  double x;
  double u_exact;

  for (int i = 0; i < N; i++) {
    x = x_min + i * dx;
    u_exact = gaussian_solution(x, t, c, g);
    I += ( * u)[i];
    R += pow(( * u)[i] - u_exact, 2.0);
    E += pow(( * u)[i], 2.0);
  }

  I = (I * dx) / (sigma * U);
  E = (E * dx) / (sigma * U * U);
  R = (R * dx) / (sigma * U * U);

  fprintf(file, "%.12f %.12f %.16f\n", I, E, R);
  return;
}

// Non uniform grid function

double g_f(double x, Grid * grid) {
  // xi = g(x) =  x - a L/(2pi) sin(2pi x / L)
  double a = grid -> a;
  double L = grid -> L;

  return x - a * L / (2 * M_PI) * sin(2 * M_PI * x / L);
}

double b(Grid * my_grid, int index) {
  // b(xi) = 1 / (1 - a cos(2pi xi / L))
  double b;
  if (my_grid -> type == UNIF) {
    b = 1.0;
  } else {
    double a = my_grid -> a;
    double L = my_grid -> L;
    double xsi = -L / 2.0 + index * my_grid -> dx;
    b = 1.0 / (1.0 - a * cos(2 * M_PI * xsi / L));
  }
  return b;
}