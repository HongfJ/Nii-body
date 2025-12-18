
#include <stdlib.h>
#include <math.h>
#include "nbody.h"

extern double const tpi;

extern double xmass[MAX_BODIES]; //mass in Earth mass, 1/333000 (3.003003003e-6) of Sun

// Angle Conversion
double zatan(double sina, double cosa);

// input: a, e, i, .....
int ele2rat(double a, double e, double ci, double w, double omega, double cm0, int ii, double xout[6]);

// cartesian to elements
int rat2ele(double xin[6], double *a, double *e, double *ci, double *w, double *omega, double *cm, int ii, int ito);


/* 轨道要素与笛卡尔坐标转换 */
int ele2rat(double a, double e, double ci, double w, double omega, double cm0, int ii, double x[6]) {
//   ii: planet ii (ii=2 for the innermost one); the last parameter=0: to the star, =1: to the barycentre of i-1 planets and the star
    double xm0 = xmass[0];
    double xm1 = xmass[ii];
    double cn = sqrt((xm0 + xm1) / pow(a, 3));
    double cm = fmod(cm0, tpi);
    if (cm <= 0.0)
    {
        cm = cm + tpi;
    }

    double ce = cm; // cm(ping jindian jiao) ; ce(pian jindian jiao)
    double ce0;
    do {
        ce0 = ce;
        ce = ce0 - (ce0 - e * sin(ce0) - cm) / (1.0 - e * cos(ce0));
    } while (fabs(ce0 - ce) > 1e-12);

    double px = cos(omega) * cos(w) - sin(omega) * sin(w) * cos(ci);
    double py = sin(omega) * cos(w) + cos(omega) * sin(w) * cos(ci);
    double pz = sin(w) * sin(ci);
    double qx = -cos(omega) * sin(w) - sin(omega) * cos(w) * cos(ci);
    double qy = -sin(omega) * sin(w) + cos(omega) * cos(w) * cos(ci);
    double qz = cos(w) * sin(ci);

    x[0] = a * (cos(ce) - e) * px + a * sqrt(1.0 - e * e) * sin(ce) * qx;
    x[1] = a * (cos(ce) - e) * py + a * sqrt(1.0 - e * e) * sin(ce) * qy;
    x[2] = a * (cos(ce) - e) * pz + a * sqrt(1.0 - e * e) * sin(ce) * qz;

    double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    x[3] = -a * a * cn / r * (sin(ce) * px - sqrt(1.0 - e * e) * cos(ce) * qx);
    x[4] = -a * a * cn / r * (sin(ce) * py - sqrt(1.0 - e * e) * cos(ce) * qy);
    x[5] = -a * a * cn / r * (sin(ce) * pz - sqrt(1.0 - e * e) * cos(ce) * qz);

    return 0;
}


int rat2ele(double x[6], double* a, double* e, double* ci, double* w, double* omega, double* cm, int ii, int ito) {
    double xm0 = xmass[0];
    double xm1 = xmass[ii];

    if (ito == 1) {
        for (int i=1; i<ii; i++){
	    xm0 =xm0 + xmass[i];
	}
    }

    double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    double eng = (x[3] * x[3] + x[4] * x[4] + x[5] * x[5]) / 2.0 - (xm0 + xm1) / r;
    *a = -(xm0 + xm1) / (2.0 * eng);

    double cn;
    cn = sqrt((xm0 + xm1) / pow(*a, 3));
    double hhx = x[1] * x[5] - x[2] * x[4];
    double hhy = -(x[0] * x[5] - x[3] * x[2]);
    double hhz = x[0] * x[4] - x[1] * x[3];
    double hh = sqrt(hhx * hhx + hhy * hhy + hhz * hhz);
    // double hu (hz in main fuction); // shengjiaodian angle ???
    // hu = hhz;
    // double hu1 (hl in main fuction); // angular momentum  ???
    // *hu1 = hh;

    double hh_squared = hh * hh;
    double tmp = hh_squared / ((xm0 + xm1) * (*a));
    if (tmp >= 1.0) {
        *e = 0.0;
    }
    else if (tmp <= 0.0) {
        *e = 1.0;
    }
    else {
        *e = sqrt(1.0 - tmp);
    }
    // 计算轨道倾角
    if (fabs(hhz / hh - 1.0) <= 1e-8) {
        *ci = 0.0;
    }
    else {
        *ci = acos(hhz / hh);
    }
    // 计算升交点赤经
    if (*ci <= 1e-8) {
        *omega = 0.0;
    }
    else {
        *omega = zatan(hhx, -hhy);
    }
    // 计算偏近点角
    double ce;
    //if ((1.0 - r / (*a)) / (*e) >= 1.0 || (*e == 0.0)) {
    if ((1.0 - r / (*a)) / (*e) >= 1.0 || (*a == r)) {
        ce = 0.0;
    }
    else if ((1.0 - r / (*a)) / (*e) <= -1.0) {
        ce = tpi / 2.0;
    }
    else {
        ce = acos((1.0 - r / (*a)) / (*e));
    }
    // 调整偏近点角方向
    if ((x[0] * x[3] + x[1] * x[4] + x[2] * x[5]) < -1e-10) {
        ce = tpi - ce;
    }
    *cm = ce - (*e) * sin(ce);
    // 计算近心点幅角
    if (*ci <= 1e-8) {
        double px = cos(ce) / r * x[0] - sin(ce) / cn / (*a) * x[3];
        double py = cos(ce) / r * x[1] - sin(ce) / cn / (*a) * x[4];
        *w = zatan(py, px);
    }
    else {
        double ee = sqrt(1.0 - (*e) * (*e));
        double pz = cos(ce) / r * x[2] - sin(ce) / cn / (*a) * x[5];
        double qz = sin(ce) / ee / r * x[2] + (cos(ce) - (*e)) / (*a) / cn / ee * x[5];
        *w = zatan(pz, qz);
    }

    return 0;
}

double zatan(double sina, double cosa) {
    double angle = atan(sina / cosa);
    if (cosa <= 0.0)               angle += tpi / 2.0;
    if (cosa >= 0.0 && sina <= 0.0) angle += tpi;
    return angle;
}


