#ifndef FD_SCHEME_H
#define FD_SCHEME_H

#include <stdlib.h>
#include "thomas.h"
#include "struct.h"


void f_E2(double** delta_u, double* ul, double const c, const int N, const double dx, const double dt) {

    for (int i=1; i<N-1; ++i){
        (*delta_u)[i] = dt *(-c / (2*dx) *(ul[i+1] - ul[i-1]));
    }
    (*delta_u)[0] = dt * (-c / (2*dx) * (ul[1] - ul[N-1]));
    (*delta_u)[N-1] = dt * (-c / (2*dx) * (ul[0] - ul[N-2]));

    return;
}

void f_ED(double** delta_u, double* ul, double const c, const int N, const double dx, const double dt) {

    for (int i=2; i<N-1; ++i){
        (*delta_u)[i] = dt * (-c) * (3*ul[i] - 6*ul[i-1] + 2*ul[i+1] + ul[i-2]) / (6*dx);
    }
    (*delta_u)[0] = dt * (-c) * (3*ul[0] - 6*ul[N-1] + 2*ul[1] + ul[N-2]) / (6*dx);
    (*delta_u)[1] = dt * (-c) * (3*ul[1] - 6*ul[0] + 2*ul[2] + ul[N-1]) / (6*dx);
    (*delta_u)[N-1] = dt * (-c) * (3*ul[N-1] - 6*ul[N-2] + 2*ul[0] + ul[N-3]) / (6*dx);

    return;
}


void f_E4(double** delta_u, double* ul, double const c, const int N, const double dx, const double dt) {
    for (int i=2; i<N-2; ++i){
        (*delta_u)[i] = dt *(-c) * (4.0/3.0 * (ul[i+1] - ul[i-1]) - 1.0/6.0 * (ul[i+2] - ul[i-2])) / (2*dx);
    }
    
    (*delta_u)[0] = dt * (-c) * (4.0/3.0 * (ul[1] - ul[N-1]) - 1.0/6.0 * (ul[2] - ul[N-2])) / (2*dx);
    (*delta_u)[1] = dt * (-c) * (4.0/3.0 * (ul[2] - ul[0]) - 1.0/6.0 * (ul[3] - ul[N-1])) / (2*dx);
    (*delta_u)[N-2] = dt * (-c) * (4.0/3.0 * (ul[N-1] - ul[N-3]) - 1.0/6.0 * (ul[0] - ul[N-4])) / (2*dx);
    (*delta_u)[N-1] = dt * (-c) * (4.0/3.0 * (ul[0] - ul[N-2]) - 1.0/6.0 * (ul[1] - ul[N-3])) / (2*dx);

    return;
}



void f_I4(double** delta_u, double* ul, double const c, const int N, const double dx, const double dt) {

    double*q = (double*) malloc(N * sizeof(double));
    
    for (int i=1; i<N-1; ++i){
        q[i] =  3.0/(4*dx) * (ul[i+1] - ul[i-1]) ;
    }
    q[0] = 3.0/(4*dx) * (ul[1] - ul[N-1]);
    q[N-1] =  3.0/(4*dx) * (ul[0] - ul[N-2]);


    solve_period_3diag(N, 1.0, 1.0/4.0, 1.0/4.0, *delta_u , q);
    
    for (int i=0; i<N; ++i){
        (*delta_u)[i] *= -c * dt;
    }
    free(q);
    return;
}

void (*select_fd_scheme(SchemeType type)) (double**, double*, double, int, double, double){
    switch (type) {
        case ED:
            return f_ED;
        case E2:
            return f_E2;
        case E4:
            return f_E4;
        case I4:
            return f_I4;
        default:
            return f_E2;
    }
}

#endif