#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nbody.h"



extern double xmass[MAX_BODIES]; //mass in Earth mass, 1/333000 (3.003003003e-6) of Sun
extern const int nbod;  // number of bodies
extern const int ndim;  // dimension of integration


/* 万有引力导数子程序 */
int f_deriv(double *x, double *y) {
    double rr3[MAX_BODIES][MAX_BODIES];
    int i, j;
    /* 速度 */
    for (i = 0; i < nbod; i++) {
        int id = ndim * i;
        y[id]   = x[id+3];
        y[id+1] = x[id+4];
        y[id+2] = x[id+5];
        y[id+3] = 0.0;
        y[id+4] = 0.0;
        y[id+5] = 0.0;
    }
    /* 计算所有两两间距的 |r|^3 */
    for (i = 0; i < nbod; i++) {
        int id = ndim*i;
        for (j = i+1; j < nbod; j++) {
            int jd = ndim*j;
            double dx = x[id] - x[jd];
            double dy = x[id+1] - x[jd+1];
            double dz = x[id+2] - x[jd+2];
            double dist_sq = dx*dx + dy*dy + dz*dz;
            rr3[i][j] = pow(dist_sq, 1.5);
            rr3[j][i] = rr3[i][j];
        }
    }
    /* 加速度 */
    for (i = 0; i < nbod; i++) {
        int id = ndim * i;
        for (j = 0; j < nbod; j++) {
            if (i == j) continue;
            int jd = ndim * j;
            y[id+3] = y[id+3] + xmass[j] * (x[jd]   - x[id])   / rr3[i][j];
            y[id+4] += xmass[j] * (x[jd+1] - x[id+1]) / rr3[i][j];
            y[id+5] += xmass[j] * (x[jd+2] - x[id+2]) / rr3[i][j];
        }
    }

    /*FILE* yhc_file = fopen("yhc.txt", "a"); 
    if (!yhc_file) { perror("yhc.txt"); return; }

    fprintf(yhc_file, "t = %g\n", t);
    for (i = 0; i < nbod; i++) {
        int id = ndim * i;
        fprintf(yhc_file, "Body %d: dx/dt = (%f, %f, %f), d2x/dt2 = (%f, %f, %f)\n",
            i, y[id], y[id + 1], y[id + 2], y[id + 3], y[id + 4], y[id + 5]);
    }

    fclose(yhc_file); */
    return 0;
}




/* 万有引力导数子程序 redundant t */
int yhc(double t, double *x, double *y) {
    double rr3[MAX_BODIES][MAX_BODIES];
    int i, j;
    /* 速度 */
    for (i = 0; i < nbod; i++) {
        int id = ndim * i;
        y[id]   = x[id+3];
        y[id+1] = x[id+4];
        y[id+2] = x[id+5];
        y[id+3] = 0.0;
        y[id+4] = 0.0;
        y[id+5] = 0.0;
    }
    /* 计算所有两两间距的 |r|^3 */
    for (i = 0; i < nbod; i++) {
        for (j = i+1; j < nbod; j++) {
            int id = ndim*i, jd = ndim*j;
            double dx = x[id] - x[jd];
            double dy = x[id+1] - x[jd+1];
            double dz = x[id+2] - x[jd+2];
            double dist_sq = dx*dx + dy*dy + dz*dz;
            rr3[i][j] = pow(dist_sq, 1.5);
            rr3[j][i] = rr3[i][j];
        }
    }
    /* 加速度 */
    for (i = 0; i < nbod; i++) {
        int id = ndim * i;
        for (j = 0; j < nbod; j++) {
            if (i == j) continue;
            int jd = ndim * j;
            y[id+3] += xmass[j] * (x[jd]   - x[id])   / rr3[i][j];
            y[id+4] += xmass[j] * (x[jd+1] - x[id+1]) / rr3[i][j];
            y[id+5] += xmass[j] * (x[jd+2] - x[id+2]) / rr3[i][j];
        }
    }

    /*FILE* yhc_file = fopen("yhc.txt", "a"); 
    if (!yhc_file) { perror("yhc.txt"); return; }

    fprintf(yhc_file, "t = %g\n", t);
    for (i = 0; i < nbod; i++) {
        int id = ndim * i;
        fprintf(yhc_file, "Body %d: dx/dt = (%f, %f, %f), d2x/dt2 = (%f, %f, %f)\n",
            i, y[id], y[id + 1], y[id + 2], y[id + 3], y[id + 4], y[id + 5]);
    }

    fclose(yhc_file); */
    return 0;
}



