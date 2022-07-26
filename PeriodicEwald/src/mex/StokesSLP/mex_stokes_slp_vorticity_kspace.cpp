#include "ewald_tools.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    if(mxGetM(prhs[0]) != 2)
        mexErrMsgTxt("psrc must be a 2xn matrix.");
    if(mxGetM(prhs[1]) != 2)
        mexErrMsgTxt("ptar must be a 2xn matrix.");
    if(mxGetM(prhs[4]) != 2)
        mexErrMsgTxt("f must be a 2xn matrix.");
    if(mxGetN(prhs[4]) != mxGetN(prhs[0]))
        mexErrMsgTxt("psrc and f must be the same size.");
    
    //Source and target points.
    double* psrc = mxGetPr(prhs[0]);
    double* ptar = mxGetPr(prhs[1]);
    int Nsrc = mxGetN(prhs[0]);    
    int Ntar = mxGetN(prhs[1]);
    
    //Ewald parameter xi
    double xi = mxGetScalar(prhs[2]);
    
    //Splitting parameter eta
    double eta = mxGetScalar(prhs[3]);
    
    //Strength vector
    double* f = mxGetPr(prhs[4]);
    
    //Number of grid intervals in each direction
    int Mx = static_cast<int>(mxGetScalar(prhs[5]));
    int My = static_cast<int>(mxGetScalar(prhs[6]));
    
    //Size of the domain
    double Lx = mxGetScalar(prhs[7]);
    double Ly = mxGetScalar(prhs[8]);
    
    //Width of the Gaussian bell curves
    double w = mxGetScalar(prhs[9]);
    
    //Number of support nodes
    int P = static_cast<int>(mxGetScalar(prhs[10]));
    
    //Grid spacing, assuming hx = hy = h
    double h = Lx/Mx;
    
    //---------------------------------------------------------------------
    //Step 1 : Spreading to the grid
    //---------------------------------------------------------------------
    //We begin by spreading the sources onto the grid. This basically means
    //that we superposition properly scaled Gaussian bells, one for
    //each source. We use fast Gaussian gridding to reduce the number of
    //exps we need to evaluate.
    
    //The function H on the grid. We need to have these as Matlab arrays
    //since we call Matlab's in-built fft2 to compute the 2D FFT.
    mxArray *fft2rhs[2],*fft2lhs[2];
    fft2rhs[0] = mxCreateDoubleMatrix(My, Mx, mxREAL);
    fft2rhs[1] = mxCreateDoubleMatrix(My, Mx, mxREAL);
    
    double* H1 = mxGetPr(fft2rhs[0]);
    double* H2 = mxGetPr(fft2rhs[1]);
    
    //This is the precomputable part of the fast Gaussian gridding.
    double* e1 = new double[P+1];
    Spread(H1, H2, e1, psrc, f, Nsrc, Lx, Ly, xi, w, eta, P, Mx, My, h);
    
    //---------------------------------------------------------------------
    //Step 2 : Frequency space filter
    //---------------------------------------------------------------------
    //We apply the appropriate k-space filter associated with the
    //Stokeslet. This involves Fast Fourier Transforms.
    
    //Call Matlab's routines to compute the 2D FFT. Other choices,
    //such as fftw could possibly be better, mainly because of the poor
    //complex data structure Matlab uses which we now have to deal with.
    mexCallMATLAB(1,&fft2lhs[0],1,&fft2rhs[0],"fft2");
    mexCallMATLAB(1,&fft2lhs[1],1,&fft2rhs[1],"fft2");

    //The output of the FFT is complex. Get pointers to the real and
    //imaginary parts of Hhat1 and Hhat2.
    
    mwSize cs = Mx*My;
    
    double* Hhat1_re = mxGetPr(fft2lhs[0]);
    double* Hhat1_im = mxGetPi(fft2lhs[0]);    
    double* Hhat2_re = mxGetPr(fft2lhs[1]);
    double* Hhat2_im = mxGetPi(fft2lhs[1]);
    
    //We cannot assume that both the real and imaginary part of the
    //Fourier transforms are non-zero.
    if(Hhat1_im == NULL) {

        Hhat1_im = (double*) mxCalloc(cs,sizeof(double));
        
        mxSetPi(fft2lhs[0],Hhat1_im);
    }
    if(Hhat2_im == NULL) {

        Hhat2_im = (double*) mxCalloc(cs,sizeof(double));
        
        mxSetPi(fft2lhs[1],Hhat2_im);
    }
   
    //Apply filter in the frequency domain. This is a completely
    //parallel operation. The FFT gives the frequency components in
    //non-sequential order, so we have to split the loops. One could
    //possibly use fast gaussian gridding here too, to get rid of the
    //exponentials, but as this loop accounts for a few percent of the
    //total runtime, this hardly seems worth the extra work.
#pragma omp parallel for
    for(int j = 0;j<Mx;j++) {
        int ptr = j*My;
        
        double k1;
        if(j <= Mx/2)
            k1 = 2.0*pi/Lx*j;
        else
            k1 = 2.0*pi/Lx*(j-Mx);
        for(int k = 0;k<=My/2;k++,ptr++) {
            
            //Hhat1 contains density component 1 convolved with Gaussians
            //Hhat2 contains density component 2 convolved with Gaussians
            double f1_re = Hhat1_re[ptr];
            double f1_im = Hhat1_im[ptr];
            double f2_re = Hhat2_re[ptr];
            double f2_im = Hhat2_im[ptr];
            
            double k2 = 2.0*pi/Ly*k;
            double Ksq = k1*k1+k2*k2;
            double e = (1.0/Ksq)*exp(-0.25*(1-eta)/(xi*xi)*Ksq);
            
            double fdotkperp_re = f1_re*k2 - f2_re*k1;
            double fdotkperp_im = f1_im*k2 - f2_im*k1;
            
            // multiplication by -i
            Hhat1_re[ptr] = fdotkperp_im*e;
            Hhat1_im[ptr] = -fdotkperp_re*e;
            
        }
        for(int k = 0;k<My/2-1;k++,ptr++) {
            
            //Hhat1 contains density component 1 convolved with Gaussians
            //Hhat2 contains density component 2 convolved with Gaussians
            double f1_re = Hhat1_re[ptr];
            double f1_im = Hhat1_im[ptr];
            double f2_re = Hhat2_re[ptr];
            double f2_im = Hhat2_im[ptr];
            
            double k2 = 2.0*pi/Ly*(k-My/2+1);
            double Ksq = k1*k1+k2*k2;
            double e = (1.0/Ksq)*exp(-0.25*(1-eta)/(xi*xi)*Ksq);
            
            double fdotkperp_re = f1_re*k2 - f2_re*k1;
            double fdotkperp_im = f1_im*k2 - f2_im*k1;
            
            // multiplication by -i
            Hhat1_re[ptr] = fdotkperp_im*e;
            Hhat1_im[ptr] = -fdotkperp_re*e;
            
        }
    }
    
    //Remove the zero frequency term.
    Hhat1_re[0] = 0;
    Hhat1_im[0] = 0;
    
    //Get rid of the old H1 and H2 arrays. They are no longer needed.
    mxDestroyArray(fft2rhs[0]);
    mxDestroyArray(fft2rhs[1]);
    
    //Do the inverse 2D FFTs. We use Matlab's inbuilt functions again.
    //Note that we could eliminate one of these IFFTs by noting that 
    //trace(grad u) = 0, but we'll keep it here for a consistency check
    mexCallMATLAB(1,&fft2rhs[0],1,&fft2lhs[0],"ifft2");
    
    //The pointer to the real part. We don't need the imaginary part.
    double* Ht1 = mxGetPr(fft2rhs[0]);
    
    if(Ht1 == NULL) {
        Ht1 = new double[Mx*My];
        memset(Ht1,0,Mx*My*sizeof(double));
    }

    //Get rid of the Hhat arrays. They are no longer needed.
    mxDestroyArray(fft2lhs[0]);
    mxDestroyArray(fft2lhs[1]);
    
    //---------------------------------------------------------------------
    //Step 3 : Evaluating the vorticity
    //---------------------------------------------------------------------
    //Here we compute the output pressure as a convolution of the Ht1
    //and Ht2 functions with a properly scaled gaussian. This is basically
    //gaussian blur, and we again use fast gaussian gridding. This
    //procedure is fully parallel.
    
    //Create the output matrix.
    plhs[0] = mxCreateDoubleMatrix(1, Ntar, mxREAL);
    double* omega = mxGetPr(plhs[0]);
    
    Gather(Ht1, 1, 1, e1, ptar, omega, Ntar, Lx, Ly, xi, w, eta, P, Mx,My, h);
        
    //Clean up
    mxDestroyArray(fft2rhs[0]);
    delete e1;
}
