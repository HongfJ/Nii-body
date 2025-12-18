#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// derivative function used in integration
int f_deriv(double *x, double *y);

//int yhc(double t, double *x, double *y) ;

/////////////////////////////////////// 
// Pruned version of the RKF7(8) algorithm for cases without time dependence.
int rkf78_no_t(double *ptr_t, double *ptr_dt, double err_tol, double *x, int n)
{
    double err_ratio = 1e-4;
    double err_now_init = 1.0;

    double *x0 = (double *)malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) {
        x0[i] = x[i];
    }
    double t0;
    t0 = *ptr_t;

    double err_now;
    err_now = err_now_init;

    while (err_now > err_tol) {
        *ptr_t = t0;
        for (int i = 0; i < n; i++) {
            x[i] = x0[i];
        }

        double *x1 = (double *)malloc(n * sizeof(double));
        double *y1 = (double *)malloc(n * sizeof(double));
        double *y2 = (double *)malloc(n * sizeof(double));
        double *y3 = (double *)malloc(n * sizeof(double));
        double *y4 = (double *)malloc(n * sizeof(double));
        double *y5 = (double *)malloc(n * sizeof(double));
        double *y6 = (double *)malloc(n * sizeof(double));
        double *y7 = (double *)malloc(n * sizeof(double));
        double *y8 = (double *)malloc(n * sizeof(double));
        double *y9 = (double *)malloc(n * sizeof(double));
        double *y10 = (double *)malloc(n * sizeof(double));
        double *y11 = (double *)malloc(n * sizeof(double));
        double *y12 = (double *)malloc(n * sizeof(double));
        double *y13 = (double *)malloc(n * sizeof(double));

        for (int i = 0; i < n; i++) {
            x1[i] = x0[i];
        }
        f_deriv(x1, y1);

        for (int i = 0; i < n; i++) {
            x1[i] = x[i] + *ptr_dt * 2.0 / 27.0 * y1[i];
        }
        f_deriv(x1, y2);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (y1[i] + 3.0 * y2[i]) / 36.0;
        }
        f_deriv(x1, y3);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (y1[i] + 3.0 * y3[i]) / 24.0;
        }
        f_deriv(x1, y4);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (y1[i] * 20.0 + (-y3[i] + y4[i]) * 75.0) / 48.0;
        }
        f_deriv(x1, y5);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (y1[i] + y4[i] * 5.0 + y5[i] * 4.0) / 20.0;
        }
        f_deriv(x1, y6);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (-y1[i] * 25.0 + y4[i] * 125.0 - y5[i] * 260.0 + y6[i] * 250.0) / 108.0;
        }
        f_deriv(x1, y7);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (y1[i] * 93.0 + y5[i] * 244.0 - y6[i] * 200.0 + y7[i] * 13.0) / 900.0;
        }
        f_deriv(x1, y8);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (y1[i] * 180.0 - y4[i] * 795.0 + y5[i] * 1408.0 - y6[i] * 1070.0 + y7[i] * 67.0 + y8[i] * 270.0) / 90.0;
        }
        f_deriv(x1, y9);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (-y1[i] * 455.0 + y4[i] * 115.0 - y5[i] * 3904.0 + y6[i] * 3110.0 - y7[i] * 171.0 + y8[i] * 1530.0 - y9[i] * 45.0) / 540.0;
        }
        f_deriv(x1, y10);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (y1[i] * 2383.0 - y4[i] * 8525.0 + y5[i] * 17984.0 - y6[i] * 15050.0 + y7[i] * 2133.0 + y8[i] * 2250.0 + y9[i] * 1125.0 + y10[i] * 1800.0) / 4100.0;
        }
        f_deriv(x1, y11);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (y1[i] * 60.0 - y6[i] * 600.0 - y7[i] * 60.0 + (y9[i] - y8[i] + 2.0 * y10[i]) * 300.0) / 4100.0;
        }
        f_deriv(x1, y12);

        for (int i = 0; i < n; i++) {
             x1[i] = x[i] + *ptr_dt * (-y1[i] * 1777.0 - y4[i] * 8525.0 + y5[i] * 17984.0 - y6[i] * 14450.0 + y7[i] * 2193.0 + y8[i] * 2550.0 + y9[i] * 825.0 + y10[i] * 1200.0 + y12[i] * 4100) / 4100.0;
        }
        f_deriv(x1, y13);

        for (int i = 0; i < n; i++) {
             x[i] = x[i] + *ptr_dt * (y6[i] * 272.0 + (y7[i] + y8[i]) * 216.0 + (y9[i] + y10[i]) * 27.0 + (y12[i] + y13[i]) * 41.0) / 840.0;
        }

        err_now = 0;
        for (int i = 0; i < n; i++) {
             double diff;
             diff = (y1[i] + y11[i] - y12[i] - y13[i]) * (*ptr_dt) * 41.0 / 810.0;
             err_now = fmax(err_now, fabs(diff));
        }

        *ptr_t = *ptr_t + *ptr_dt;

        if (err_now > err_tol) {
             *ptr_dt = *ptr_dt / 2.0;
        }

        if (err_now < err_tol * err_ratio) {
             *ptr_dt = *ptr_dt * 2.0;
        }

        free(x1);
	x1 = NULL;
        free(y1);
	y1 = NULL;
        free(y2);
	y2 = NULL;
        free(y3);
	y3 = NULL;
        free(y4);
	y4 = NULL;
        free(y5);
	y5 = NULL;
        free(y6);
	y6 = NULL;
        free(y7);
	y7 = NULL;
        free(y8);
	y8 = NULL;
        free(y9);
	y9 = NULL;
        free(y10);
	y10 = NULL;
        free(y11);
	y11 = NULL;
        free(y12);
	y12 = NULL;
        free(y13);
	y13 = NULL;
    }

    free(x0);
    x0 = NULL;

    return 0;
}




/////////////////////////////////////////////
// example derivative function used in integration
int yhc(double t, double *x, double *y);
//
//////////////////// Full RKF7(8) algorithm
int rkf78(double *ptr_t, double *ptr_dt, double err_tol, double *x, int n) {
    // Runge Kutta Fehlberg 7(8) numerical integration method
    
    /*
     t  (*ptr_t):  start time of this integration (t_init)
     dt (*ptr_dt): initial step size attempted to use (may fail due to large err)
     *
     ptr_t:        pointer to t
     ptr_dt:       pointer to dt
     err_tol:      truncation error tolerance [non-dimensional]
     *x:           integration vector at time = ti
     n:            number of differential equations (dimension of x)
     *
     RETURN:  *x with updated values at time tf.
     *
     */

    // local constants
    // 
    // criteria in the tiny err case to increase step size
    double err_ratio = 1e-4; 
    // set manually err_now to big value to start the loop
    double err_now_init = 1.0;


    //////////////////////////////
    // Create a copy of initial x and t, now x0 and t0
    //
    // allocate arrays
    double *x0 = (double *)malloc(n * sizeof(double));
    // copy initial values
    for (int i = 0; i < n; i++) 
    {
            x0[i] = x[i];
    }
    //
    // integration start time at each step
    double t0;
    t0 = *ptr_t;


    ////////////////////////////
    //  
    // err for the working step size
    double err_now;
    // set it to 1.0 at the beginning
    err_now = err_now_init;  


    /////////////////////////////////////////////////
    // Main loop
    //
    while (err_now > err_tol) 
    {

        ///////////////////
	// re-initialize the t, t1, *x
	*ptr_t = t0;
        double t1;
	t1 = *ptr_t;
	//
        for (int i = 0; i < n; i++) 
	{
            x[i] = x0[i];
        }
        
	// allocate arrays of x and integration coefficients
        double *x1 = (double *)malloc(n * sizeof(double));
        double *y1 = (double *)malloc(n * sizeof(double));
        double *y2 = (double *)malloc(n * sizeof(double));
        double *y3 = (double *)malloc(n * sizeof(double));
        double *y4 = (double *)malloc(n * sizeof(double));
        double *y5 = (double *)malloc(n * sizeof(double));
        double *y6 = (double *)malloc(n * sizeof(double));
        double *y7 = (double *)malloc(n * sizeof(double));
        double *y8 = (double *)malloc(n * sizeof(double));
        double *y9 = (double *)malloc(n * sizeof(double));
        double *y10 = (double *)malloc(n * sizeof(double));
        double *y11 = (double *)malloc(n * sizeof(double));
        double *y12 = (double *)malloc(n * sizeof(double));
        double *y13 = (double *)malloc(n * sizeof(double));


	//////////////////
	// set x1 = x0
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x0[i];
        }
        // calc y1 use derivative function
        yhc(t1, x1, y1);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * 2.0 / 27.0 * y1[i];
        }
        t1 = *ptr_t + *ptr_dt * 2.0 / 27.0;
	//
        yhc(t1, x1, y2);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (y1[i] + 3.0 * y2[i]) / 36.0;
        }
        t1 = *ptr_t + *ptr_dt / 9.0;
	//
        yhc(t1, x1, y3);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (y1[i] + 3.0 * y3[i]) / 24.0;
        }
        t1 = *ptr_t + *ptr_dt / 6.0;
	//
        yhc(t1, x1, y4);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (y1[i] * 20.0 + (-y3[i] + y4[i]) * 75.0) / 48.0;
        }
        t1 = *ptr_t + *ptr_dt * 5.0 / 12.0;
	//
        yhc(t1, x1, y5);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (y1[i] + y4[i] * 5.0 + y5[i] * 4.0) / 20.0;
        }
        t1 = *ptr_t + *ptr_dt * 1.0 / 2.0;
	//
        yhc(t1, x1, y6);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (-y1[i] * 25.0 + y4[i] * 125.0 - y5[i] * 260.0 + y6[i] * 250.0) / 108.0;
        }
        t1 = *ptr_t + *ptr_dt * 5.0 / 6.0;
	//
        yhc(t1, x1, y7);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (y1[i] * 93.0 + y5[i] * 244.0 - y6[i] * 200.0 + y7[i] * 13.0) / 900.0;
        }
        t1 = *ptr_t + *ptr_dt * 1.0 / 6.0;
	//
        yhc(t1, x1, y8);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (y1[i] * 180.0 - y4[i] * 795.0 + y5[i] * 1408.0 - y6[i] * 1070.0 + y7[i] * 67.0 + y8[i] * 270.0) / 90.0;
        }
        t1 = *ptr_t + *ptr_dt * 2.0 / 3.0;
	//
        yhc(t1, x1, y9);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (-y1[i] * 455.0 + y4[i] * 115.0 - y5[i] * 3904.0 + y6[i] * 3110.0 - y7[i] * 171.0 + y8[i] * 1530.0 - y9[i] * 45.0) / 540.0;
        }
        t1 = *ptr_t + *ptr_dt * 1.0 / 3.0;
	//
        yhc(t1, x1, y10);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (y1[i] * 2383.0 - y4[i] * 8525.0 + y5[i] * 17984.0 - y6[i] * 15050.0 + y7[i] * 2133.0 + y8[i] * 2250.0 + y9[i] * 1125.0 + y10[i] * 1800.0) / 4100.0;
        }
        t1 = *ptr_t + *ptr_dt;
	//
        yhc(t1, x1, y11);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (y1[i] * 60.0 - y6[i] * 600.0 - y7[i] * 60.0 + (y9[i] - y8[i] + 2.0 * y10[i]) * 300.0) / 4100.0;
        }
        t1 = *ptr_t;
	//
        yhc(t1, x1, y12);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x1[i] = x[i] + *ptr_dt * (-y1[i] * 1777.0 - y4[i] * 8525.0 + y5[i] * 17984.0 - y6[i] * 14450.0 + y7[i] * 2193.0 + y8[i] * 2550.0 + y9[i] * 825.0 + y10[i] * 1200.0 + y12[i] * 4100) / 4100.0;
        }
        t1 = *ptr_t + *ptr_dt;
	//
        yhc(t1, x1, y13);

        //////////
        for (int i = 0; i < n; i++) 
	{
            x[i] = x[i] + *ptr_dt * (y6[i] * 272.0 + (y7[i] + y8[i]) * 216.0 + (y9[i] + y10[i]) * 27.0 + (y12[i] + y13[i]) * 41.0) / 840.0;
        }

        /////////
	// find the largest error
        err_now = 0;
        //
        for (int i = 0; i < n; i++)
	{
            double diff;
            diff = (y1[i] + y11[i] - y12[i] - y13[i]) * (*ptr_dt) * 41.0 / 810.0;
            err_now = fmax(err_now, fabs(diff));
        }

        *ptr_t = *ptr_t + *ptr_dt;

        if (err_now > err_tol) 
	{
            *ptr_dt = *ptr_dt / 2.0;
        }

        if (err_now < err_tol * err_ratio)
	{
            *ptr_dt = *ptr_dt * 2.0;
        }
	
	//
        /////////////////////
        // free the allocated arrays
        free(x1);
	x1 = NULL;
        free(y1);
	y1 = NULL;
        free(y2);
	y2 = NULL;
        free(y3);
	y3 = NULL;
        free(y4);
	y4 = NULL;
        free(y5);
	y5 = NULL;
        free(y6);
	y6 = NULL;
        free(y7);
	y7 = NULL;
        free(y8);
	y8 = NULL;
        free(y9);
	y9 = NULL;
        free(y10);
	y10 = NULL;
        free(y11);
	y11 = NULL;
        free(y12);
	y12 = NULL;
        free(y13);
	y13 = NULL;
    }

    free(x0);
    x0 = NULL;

    return 0;
}




