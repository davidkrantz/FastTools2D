#include "mex.h"
#include "omp.h"
#include "math.h"
#define pi 3.1415926535897932385

void fgg_exp(double*,double,int,double,double,double,int);

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* Called from matlab as
     * H = mex_se_K0_2p_quadrature(ptar,Htilde,Z0,h,L,afac,M,P)
     *
     * Input:
     *      ptar:   target points in an 2xn array
     *      Htilde: matrix of size MxM
     *      Z0,h,L,afac,M,P: ewald parameters
     *
     * Output:
     *      u: the approximated integrals of size 1xn
     */
    
    /* Check nbr of inputs and outputs */
    if(nrhs != 10)
        mexErrMsgIdAndTxt("mex_se_K0_2p_quadrature:invalidNumInputs",
                "Eight inputs required.");
    if(nlhs != 1)
        mexErrMsgIdAndTxt("mex_se_K0_2p_quadrature:maxlhs",
                "Too many output arguments.");
    
    /* Declare input variables */
    double *ptar = mxGetPr(prhs[0]);              // Targets in a 2xNtar array
    int Ntar = mxGetN(prhs[0]);                   // Nbr of targets
    double *Htilde = mxGetPr(prhs[1]);            // Matrix of size MxM
    double *Z0 = mxGetPr(prhs[2]);                // Pre-computed matrix of size PxP for Fast Gaussian Gridding
    double h = mxGetScalar(prhs[3]);              // Grid spacing uniform grid
    double Lx = mxGetScalar(prhs[4]);              // Length of one side of reference cell
    double Ly = mxGetScalar(prhs[5]);              // Length of one side of reference cell
    double afac = mxGetScalar(prhs[6]);           // alpha-factor of window function width
    int Mx = mxGetScalar(prhs[7]);                 // Number of uniform disc. points in x direction
    int My = mxGetScalar(prhs[8]);                 // Number of uniform disc. points in y direction
    int P = mxGetScalar(prhs[9]);                 // Nbr of points in window func. support
    /*
     * Note that elements in Z0 are accessed by Z0[j*P+k] for column j, row k.
     * Elements of Htilde are accessed by Htilde[j*M+k] for column j, row k
     */
    
    /* Check that psrc, fsrc and Z0 have the correct size */
    if(mxGetM(prhs[0]) != 2)
        mexErrMsgTxt("mex_se_K0_2p_quadrature:psrc must be a 2xn matrix.");
    if(mxGetM(prhs[1]) != My || mxGetN(prhs[1]) != Mx)
        mexErrMsgTxt("mex_se_K0_2p_quadrature:Htilde must be an (My)x(Mx) matrix.");
    if(mxGetM(prhs[2]) != P || mxGetN(prhs[2]) != P)
        mexErrMsgTxt("mex_se_K0_2p_quadrature:Z0 must be a PxP matrix.");
    
    /* Declare output variables */
    plhs[0] = mxCreateDoubleMatrix(1, Ntar, mxREAL);
    double* u = mxGetPr(plhs[0]);
    
    /* Compute trapezoidal weights. Equal since periodic */
    double wmat = h*h;
    
    /* Go through all sources */
#pragma omp parallel for
    for(int m = 0; m < Ntar; m++) {
        
        double xm = ptar[2*m];
        double ym = ptar[2*m+1];
        
        /* Find closest grid point */
        int jm = floor((xm+Lx/2.0)/h);
        int km = floor((ym+Ly/2.0)/h);
        
        /* Find offset from closest grid point to (xn,yn) */
        double txm = xm - (jm*h-Lx/2.0);
        double tym = ym - (km*h-Ly/2.0);
        
        double Zm = exp(-afac*(txm*txm+tym*tym));
        double* Zx = malloc(P*sizeof(double));
        fgg_exp(Zx,xm,jm,afac,Lx,h,P);
        double* Zy = malloc(P*sizeof(double));
        fgg_exp(Zy,ym,km,afac,Ly,h,P);
        
        double am = afac/pi;
        
        for (int j=0; j<P; j++) {
            /* Find global index for x-direction */
            int idxj = jm + j - (P/2-1);
            /* Wrap around for periodicity */
            if (idxj < 0)
                idxj += Mx;
            
            idxj = idxj % Mx;
            
            
            for (int k=0; k<P; k++) {
                /* Find global index for y-direction */
                int idxk = km + k - (P/2-1);
                /* Wrap around for periodicity */
                if (idxk < 0)
                    idxk += My;
                
                idxk = idxk % My;
                
                
                /* Compute u: idxj, y:idxk */
                u[m] += Htilde[idxj*My+idxk]*wmat*am*Z0[j*P+k]*Zm*Zx[j]*Zy[k];
            }
        }
        
        free(Zx);
        free(Zy);
        
    }
}

void fgg_exp(double* Z, double x, int j, double af, double L, double h, int P) {
    /* Computes the Zx and Zy contributions for the Fast Gaussian Gridding.
     * Z has size Px1 */
    
    double t = x-(h*j-L/2.0); //Note that this assumes a square box...
    
    double q = exp(2*af*h*t);
    Z[0] = pow(q,-P/2.0+1);
    for (int k=1; k<P; k++) {
        Z[k] = q*Z[k-1];
    }
}