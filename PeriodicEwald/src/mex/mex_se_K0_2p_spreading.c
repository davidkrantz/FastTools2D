#include "mex.h"
#include "omp.h"
#include "math.h"
#define pi 3.1415926535897932385

void fgg_exp(double*,double,int,double,double,double,int);

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* Called from matlab as
     * H = mex_se_K0_2p_spreading(psrc,fsrc,ewaldparam)
     *
     * Input:
     *      psrc: a 2xn array for n source points (x and y coord.)
     *      fsrc: a 1xn array for n source forces
     *      ewaldparam: a struct containing all spectral Ewald parameters
     * Output:
     *      H: a MxM matrix
     */
    
    /* Check nbr of inputs and outputs */
    if(nrhs != 10)
        mexErrMsgIdAndTxt("mex_se_K0_2p_spreading:invalidNumInputs",
                "Eight inputs required.");
    if(nlhs != 1)
        mexErrMsgIdAndTxt("mex_se_K0_2p_spreading:maxlhs",
                "Too many output arguments.");
    
    /* Declare input variables */
    double *psrc = mxGetPr(prhs[0]);              // Sources in a 2xNsrc array
    double *fsrc = mxGetPr(prhs[1]);              // Forces at source points in a 1xNsrc array
    int Nsrc = mxGetN(prhs[0]);                   // Nbr of sources
    int Mx = mxGetScalar(prhs[2]);                 // Number of uniform disc. points in x direction
    int My = mxGetScalar(prhs[3]);                 // Number of uniform disc. points in y direction
    int P = mxGetScalar(prhs[4]);                 // Nbr of points in window func. support
    double h = mxGetScalar(prhs[5]);              // Grid spacing uniform grid
    double afac = mxGetScalar(prhs[6]);           // alpha-factor of window function width
    double Lx = mxGetScalar(prhs[7]);              // Length of x side of reference cell
    double Ly = mxGetScalar(prhs[8]);              // Length of y side of reference cell
    double *Z0 = mxGetPr(prhs[9]);                // Pre-computed matrix of size PxP for Fast Gaussian Gridding
// Note that elements in Z0 are accessed by Z0[j*P+k] for column j, row k.
    
    /* Check that psrc, fsrc and Z0 have the correct size */
    if(mxGetM(prhs[0]) != 2)
        mexErrMsgTxt("mex_se_K0_2p_spreading:psrc must be a 2xn matrix.");
    if(mxGetM(prhs[1]) != 1)
        mexErrMsgTxt("mex_se_K0_2p_spreading:fsrc must be a 1xn matrix.");
    if(mxGetM(prhs[9]) != P || mxGetN(prhs[9]) != P)
        mexErrMsgTxt("mex_se_K0_2p_spreading:Z0 must be a PxP matrix.");
    
    /* Declare output variables */
    plhs[0] = mxCreateDoubleMatrix(My, Mx, mxREAL);
    double* H = mxGetPr(plhs[0]);
    
    /*
     * Note that we can't use openMP easily here, since the indices wrap around.
     * We could perhaps lock parts of the matrix, but this is not implemented...
     */
    
    /* Go through all sources */
    for(int n = 0; n<Nsrc; n++) {
        
        double xn = psrc[2*n];
        double yn = psrc[2*n+1];
        double fn = fsrc[n];
        
        /* Find closest grid point */
        int jn = floor((xn+Lx/2.0)/h);
        int kn = floor((yn+Ly/2.0)/h);
        
        /* Find offset from closest grid point to (xn,yn) */
        double txn = xn - (jn*h-Lx/2.0);
        double tyn = yn - (kn*h-Ly/2.0);
        
        double Zn = exp(-afac*(txn*txn+tyn*tyn));
        double* Zx = malloc(P*sizeof(double));
        fgg_exp(Zx,xn,jn,afac,Lx,h,P);
        double* Zy = malloc(P*sizeof(double));
        fgg_exp(Zy,yn,kn,afac,Ly,h,P);
        
        double an = afac*fn/pi;
        
        for (int j=0; j<P; j++) {
            /* Find global index for x-direction */
            int idxj = jn + j - (P/2-1);
            /* Wrap around for periodicity */
            if (idxj < 0)
                idxj += Mx;
            
            idxj = idxj % Mx;
            
            for (int k=0; k<P; k++) {
                /* Find global index for y-direction */
                int idxk = kn + k - (P/2-1);
                /* Wrap around for periodicity */
                if (idxk < 0)
                    idxk += My;
                
                idxk = idxk % My;
                
                
                /* Compute H contribution x: idxj, y:idxk */
                H[idxj*My+idxk] += an*Z0[j*P+k]*Zn*Zx[j]*Zy[k];
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