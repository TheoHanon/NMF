#ifndef STRUCT_H
#define STRUCT_H

#include <stdio.h>

#include <math.h>

#include <stdlib.h>

#include <string.h>

typedef enum {
  UNIF,
  NONUNIF
}
GridType;
typedef enum {
  ED,
  E2,
  E4,
  I4
}
SchemeType;
typedef enum {
  GAUSSIAN,
  WAVEPACKET
}
InitialType;

GridType get_grid_type(const char * type);
SchemeType get_scheme_type(const char * type);
InitialType get_initial_type(const char * type);

typedef struct Parameter {
  int N;
  double L;
  double dx;
  double a;
  double dt;
  double beta;
  double c;
  int N_step;
  double U;
  GridType type;
  SchemeType scheme_type;
}
Parameter;

typedef struct Grid {
  double dx;
  int N;
  GridType type;
  double * b;
  double a;
  double L;
}
Grid;

typedef struct Gaussian {
  double U;
  double mu;
  double sigma;
  double xmin;
  double xmax;
  double L;
}
Gaussian;

typedef struct Scheme {
  double dt;
  double dx;
  int N;
  SchemeType type;
}
Scheme;

Gaussian * init_gaussian(double U, double mu, double dx, double L);
void initGaussianArray(double ** arr, int n, Grid * grid, struct Gaussian * g);
double gaussian_solution(const double x,
  const double t,
    const double c, struct Gaussian * g);
void initWavePacketArray(double ** arr, int n, int p, Gaussian * g, Grid * grid);
double gaussian(double x, struct Gaussian * g);

void save(FILE * fptr,
  const double * u,
    const double t,
      const int N);

double g_f(double x, Grid * grid);
double b(Grid * my_grid, int index);

void diagnose(FILE * file, double ** u, double t,
  const int N,
    const double dx, struct Gaussian * g,
      const double c);

Grid * init_grid(const double dx,
  const int N,
    const GridType type,
      const double a,
        const double L);
void free_grid(Grid * g);
void save_grid(Grid * grid, char * file_name);

Parameter * init_parameter(int N, double L, double dx, double a, double dt, double beta, double c, int N_step, double U, GridType type, SchemeType scheme_type);
void free_parameter(Parameter * param);
void save_parameter(Parameter * param, char * file_name);

Scheme * init_scheme(const double dt,
  const double dx,
    const int N, SchemeType scheme_type);
void free_scheme(Scheme * s);

#endif