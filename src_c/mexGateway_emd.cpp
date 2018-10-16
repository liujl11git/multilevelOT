// This is a MATLAB interface to call the function in "emdPDHG.c".
// This file is edited by Jialin Liu @ UCLA (liujl11@math.ucla.edu).

#include <math.h>
#include <stdlib.h>
#include "mex.h"
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdio.h>


void emdPDHG(int p, double *density, double tau, double sigma, double error, int maxIters, int n1, int n2, double *ux, double *uy,double *px,double *py, double *psix, double *psiy, double *kpotential, double * iters, double * ffttime, double * calctime);


/* The gateway function */
// // [ux,uy,px,py,psix,psiy,potential,itrs,ffttime,calctime] = emdPDHG(density, tau, sigma, error, maxIters, initialmx, initialmy,initialpx,initialpy, p)
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

// inputs
int n1 = mxGetM(prhs[0]);
int n2 = mxGetN(prhs[0]);
double *density = mxGetPr(prhs[0]);
double tau = mxGetScalar(prhs[1]);
double sigma = mxGetScalar(prhs[2]);
double error = mxGetScalar(prhs[3]);
int maxIters = (int)mxGetScalar(prhs[4]);

/* code here */
double *ux = mxGetPr(prhs[5]);
double *uy = mxGetPr(prhs[6]);
double *px = mxGetPr(prhs[7]);
double *py = mxGetPr(prhs[8]);
int p = (int)mxGetScalar(prhs[9]);
double * psix = (double *)calloc(n1*n2,sizeof(double));
double * psiy = (double *)calloc(n1*n2,sizeof(double));
double *kpotential=(double *)calloc(n1*n2,sizeof(double));
double * iters = (double *)calloc(1,sizeof(double));
double * ffttime = (double *)calloc(1,sizeof(double));
double * calctime = (double *)calloc(1,sizeof(double));
emdPDHG(p, density, tau, sigma, error, maxIters, 
    n1, n2, ux, uy, px, py, psix, psiy,kpotential,iters,ffttime,calctime);

//outputs
plhs[0] = mxCreateDoubleMatrix(n1, n2, mxREAL);
plhs[1] = mxCreateDoubleMatrix(n1, n2, mxREAL);
memcpy(mxGetPr(plhs[0]), ux, sizeof(double)*n1*n2);
memcpy(mxGetPr(plhs[1]), uy, sizeof(double)*n1*n2);
plhs[2] = mxCreateDoubleMatrix(n1, n2, mxREAL);
plhs[3] = mxCreateDoubleMatrix(n1, n2, mxREAL);
memcpy(mxGetPr(plhs[2]), px, sizeof(double)*n1*n2);
memcpy(mxGetPr(plhs[3]), py, sizeof(double)*n1*n2);
plhs[4] = mxCreateDoubleMatrix(n1, n2, mxREAL);
plhs[5] = mxCreateDoubleMatrix(n1, n2, mxREAL);
memcpy(mxGetPr(plhs[4]), psix, sizeof(double)*n1*n2);
memcpy(mxGetPr(plhs[5]), psiy, sizeof(double)*n1*n2);
plhs[6] = mxCreateDoubleMatrix(n1, n2, mxREAL);
memcpy(mxGetPr(plhs[6]), kpotential, sizeof(double)*n1*n2);
plhs[7] = mxCreateDoubleMatrix(1, 1, mxREAL);
memcpy(mxGetPr(plhs[7]), iters, sizeof(double));
plhs[8] = mxCreateDoubleMatrix(1, 1, mxREAL);
memcpy(mxGetPr(plhs[8]), ffttime, sizeof(double));
plhs[9] = mxCreateDoubleMatrix(1, 1, mxREAL);
memcpy(mxGetPr(plhs[9]), calctime, sizeof(double));

free(psix);
free(psiy);
free(kpotential);
free(iters);
free(calctime);
}


