#ifndef RK4_H
#define RK4_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>


#include "fd_scheme.h"
#include "struct.h"


double alpha[] = {0.0, 1.0/2.0, 1.0/2.0, 1.0};
double gamma[] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};


void RK4(double *us, double** p_ul, double** p_u, double** p_delta_u,
                     double* p_t, const double c, struct Scheme* s, Grid* grid){

    // delta_u must be full of zeros !
    // us = u !!!

    double * delta_u = *p_delta_u;
    double * ul = *p_ul;
    double * u = *p_u; 

    // initialize a function pointer 
    void (*f)(double**, double*, double, int, double, double);
    f = select_fd_scheme(s->type);
    
    for (int k = 0; k < 4; k++){
        for (int i=0; i < s->N; i++) {
            ul[i] = b(grid, i) * (us[i] + alpha[k]*delta_u[i]);
        }

        (*f)(&delta_u, ul, c, s->N, s->dx, s->dt);

        for (int i=0; i < s->N; i++) {
            u[i] = u[i] + gamma[k]*delta_u[i];
        }
        *p_t = *p_t + gamma[k]*s->dt;
    }

    return;

}

void convection_solve(char* data_file_name, char* diagnose_file_name, double** p_u0, const double N_step, const double c, Grid* grid, struct Scheme* s, struct Gaussian* g) {

    FILE* data_pf = fopen(data_file_name, "w"); 
    FILE* diagnose_pf = fopen(diagnose_file_name, "w");

    double* u = *p_u0;
    double* ul = (double*) calloc(s->N, sizeof(double));
    double* us = (double*) malloc(s -> N * sizeof(double));
    double* delta_u = (double*) malloc(s -> N * sizeof(double));
    double t = 0.0;

    diagnose(diagnose_pf, &u, t, s -> N, s -> dx, g,c);
    save(data_pf, u, t, s -> N);
    for (int i=0; i < N_step; ++i){
        memcpy(us, u, s -> N*sizeof(double));
        memset(delta_u, 0.0, s -> N*sizeof(double));
        RK4(us, &ul, &u, &delta_u, &t, c, s, grid);
        diagnose(diagnose_pf, &u, t, s -> N, grid -> dx, g, c);
        save(data_pf, u, t, s ->N);
    }
    
 
    
    fclose(data_pf);
    fclose(diagnose_pf);
    free(ul);
    free(us);
    free(delta_u);

    return;

}


#endif