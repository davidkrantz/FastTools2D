#include "mex.h"
#include <bessellib.h>
#include "omp.h"
#include <gsl/gsl_sf_expint.h>


/* This function calls the C function K0xy from bessellib.c */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
/* Called from matlab as K0 = compK0(z) */
    
/* Check inputs and outputs */
if(nrhs != 2) {
    mexErrMsgIdAndTxt("compK0:nrhs",
                      "Two input required.");
}
if(nlhs != 1) {
    mexErrMsgIdAndTxt("compK0:nlhs",
                      "One output required.");
}

/* Declare variables */
double* x = mxGetPr(prhs[0]);                                              // The indata x
int ndatax = mxGetM(prhs[0]);                                              // Size of indata x
double* y = mxGetPr(prhs[1]);                                              // The indata y
int ndatay = mxGetM(prhs[1]);                                              // Size of indata y

plhs[0] = mxCreateDoubleMatrix(ndatax*ndatay,1,mxREAL);                    // Create output data array
double* K0 = mxGetPr(plhs[0]);                                           // Point to output data by K0z

// For each target in z, call besselk0 from mathint.c
#pragma omp parallel for
for (int i=0; i<ndatax; i++) {
    for (int j=0; j<ndatay; j++) {
//         K0[i + j*ndatax] = gsl_sf_expint_E1(x[i]);
        
        int alg = 0;
        int jmax = choose_algorithm(x[i], y[j], &alg);
        
//         int alg = 1;
//         int jmax = ceil(15.-x[i]);
//         if (jmax < 1) jmax = 2;
        
        K0[i + j*ndatax] = K0xy(x[i],y[j],alg,jmax);
        
//         if (alg != 2) {
//             mexPrintf("i=%d, j=%d, alg = %d, jmax = %d\n",i,j,alg, jmax);
//         }
    }
}


}