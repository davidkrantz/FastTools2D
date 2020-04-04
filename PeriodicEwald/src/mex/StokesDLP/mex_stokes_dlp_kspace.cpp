#include "mex.h"
#include <math.h>
#include <omp.h>
#include <string.h>
#define pi 3.1415926535897932385

           
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    if(nrhs != 12)
        mexErrMsgTxt("Incorrect number of input parameters");
    
    if(mxGetM(prhs[0]) != 2)
        mexErrMsgTxt("psrc must be a 2xn matrix.");
    if(mxGetM(prhs[1]) != 2)
        mexErrMsgTxt("ptar must be a 2xn matrix.");
    if(mxGetM(prhs[4]) != 2)
        mexErrMsgTxt("f must be a 2xn matrix.");
    if(mxGetM(prhs[5]) != 2)
        mexErrMsgTxt("n must be a 2xn matrix.");
    if(mxGetN(prhs[4]) != mxGetN(prhs[0]))
        mexErrMsgTxt("psrc and f must be the same size.");
    if(mxGetN(prhs[5]) != mxGetN(prhs[0]))
        mexErrMsgTxt("psrc and n must be the same size.");
    
    //The points.
    double* psrc = mxGetPr(prhs[0]);  
    int Nsrc = mxGetN(prhs[0]);
    double* ptar = mxGetPr(prhs[1]);  
    int Ntar = mxGetN(prhs[1]);
    
    
    //The Ewald parameter xi
    double xi = mxGetScalar(prhs[2]);
    //The splitting parameter eta
    double eta = mxGetScalar(prhs[3]);
    //The Stresslet vectors.
    double* f = mxGetPr(prhs[4]);
    double* n = mxGetPr(prhs[5]);
    
    int M_x = static_cast<int>(mxGetScalar(prhs[6]));
    int M_y = static_cast<int>(mxGetScalar(prhs[7]));
    
    //The length of the domain
    double L_x = mxGetScalar(prhs[8]);
    double L_y = mxGetScalar(prhs[9]);
    
    //The width of the Gaussian bell curves. (unnecessary?)
    double w = mxGetScalar(prhs[10]);
    
    //The width of the Gaussian bell curves on the grid.
    int P = static_cast<int>(mxGetScalar(prhs[11]));
    
    //The grid spacing
    double h_x = L_x/M_x;
    double h_y = L_y/M_y;
    
    //---------------------------------------------------------------------
    //Step 1 : Spreading to the grid
    //---------------------------------------------------------------------
    //We begin by spreading the sources onto the grid. This basically means
    //that we superposition properly scaled Gaussian bells, one for
    //each source. We use fast Gaussian gridding to reduce the number of
    //exps we need to evaluate.
  
    //The function H on the grid. We need to have these as Matlab arrays 
    //since we call Matlab's in-built fft2 to compute the 2D FFT.
    mxArray *fft2rhs[4],*fft2lhs[4];
    fft2rhs[0] = mxCreateDoubleMatrix(M_y, M_x, mxREAL);
    double* H1 = mxGetPr(fft2rhs[0]);
    fft2rhs[1] = mxCreateDoubleMatrix(M_y, M_x, mxREAL);
    double* H2 = mxGetPr(fft2rhs[1]);
    fft2rhs[2] = mxCreateDoubleMatrix(M_y, M_x, mxREAL);
    double* H3 = mxGetPr(fft2rhs[2]);
    fft2rhs[3] = mxCreateDoubleMatrix(M_y, M_x, mxREAL);
    double* H4 = mxGetPr(fft2rhs[3]);
    
    //This is the precomputable part of the fast Gaussian gridding.
    double* e1 = new double[P+1];
    double tmp = -2*xi*xi/eta*h_x*h_y;
    for(int j = -P/2;j<=P/2;j++)
        e1[j+P/2] = exp(tmp*j*j);
    
    double TOL = 1e-13;
    
    //Spreading the sources to the grid is not a completely parallel 
    //operation. We use the simple approach of locking the column of the 
    //matrix we are working on currently. This might not be optimal but it 
    //is simple to implement.
    omp_lock_t* locks = new omp_lock_t[M_x];
    for(int j = 0;j<M_x;j++)
        omp_init_lock(&locks[j]);
    
    //We use OpenMP for simple parallelization.
    #pragma omp parallel for 
    for(int k = 0;k<Nsrc;k++) {

        //The Gaussian bells are translation invariant. We exploit this 
        //fact to avoid blow-up of the terms and the numerical instability
        //that follows. (px,py) is the center of the bell with the original
        //grid-alignment but close to the origin.
        
        double px = psrc[2*k] - h_x*floor(psrc[2*k]/h_x);
        double py = psrc[2*k+1] - h_y*floor(psrc[2*k+1]/h_y);

        double Lhalf_x = static_cast<double>(L_x/2.0);
        double Lhalf_y = static_cast<double>(L_y/2.0);
        
        //(mx,my) is the center grid point in H1 and H2.        
        int mx;
        double tmpp = static_cast<double>((psrc[2*k]+Lhalf_x)/h_x);
        int tmppint = static_cast<int>(floor(tmpp));
        int tmppint2 = static_cast<int>(round(tmpp));
        int tmppint3 = static_cast<int>(ceil(tmpp-1));
        if(fabs(tmpp-tmppint2)<TOL) {
            tmppint = tmppint2;
            mx = static_cast<int>(tmppint-P/2);
        }else
        {
            if(fabs(px)>TOL)
                mx = static_cast<int>(tmppint3-P/2);
            else
                mx = static_cast<int>(tmppint-P/2);
        }
        
        int my;
        double tmppy = static_cast<double>((psrc[2*k+1]+Lhalf_y)/h_y);
        int tmppinty = static_cast<int>(floor(tmppy));
        int tmppint2y = static_cast<int>(round(tmppy));
        int tmppint3y = static_cast<int>(ceil(tmppy-1));
        if(fabs(tmppy-tmppint2y)<TOL) {
            tmppinty = tmppint2y;
            my = static_cast<int>(tmppinty-P/2);            
        }
        else
        {
            if(fabs(py)>TOL)
                my = static_cast<int>(tmppint3y-P/2);
            else
                my = static_cast<int>(tmppinty-P/2);
        }
        
        if (fabs(px - h_x) < TOL)
            px = 0;
        
        if (fabs(py - h_y) < TOL)
            py = 0;
        
        //Some auxillary quantities for the fast Gaussian gridding.
        double tmp = -2*xi*xi/eta;
        double ex = exp(tmp*(px*px+py*py + 2*w*px));
        double e4y = exp(2*tmp*w*py);
        double e3x = exp(-2*tmp*h_x*px);
        double e3y = exp(-2*tmp*h_y*py);
        
        //We add the Gaussians column by column, and lock the one we are
        //working on to avoid race conditions.
        
        for(int x = 0;x<P+1;x++) {
            double ey = ex*e4y*e1[x];
            int xidx = ((x+mx+M_x)%M_x)*M_y;
            
            omp_set_lock(&locks[(x+mx+M_x)%M_x]);
            if(my >= 0 && my < M_y-P-1) {
                int idx = my+xidx;
                
                for(int y = 0;y<P+1;y++,idx++) {
                    double tmp = ey*e1[y];
                    H1[idx] += tmp*n[2*k]*f[2*k];
                    H2[idx] += tmp*n[2*k+1]*f[2*k];
                    H3[idx] += tmp*n[2*k]*f[2*k+1];
                    H4[idx] += tmp*n[2*k+1]*f[2*k+1];
                    ey *= e3y;
                }     
            }else{
                for(int y = 0;y<P+1;y++) {
                    double tmp = ey*e1[y];
                    int idx = ((y+my+M_y)%M_y)+xidx;
                    H1[idx] += tmp*n[2*k]*f[2*k];
                    H2[idx] += tmp*n[2*k+1]*f[2*k];
                    H3[idx] += tmp*n[2*k]*f[2*k+1];
                    H4[idx] += tmp*n[2*k+1]*f[2*k+1];
                    ey *= e3y;
                }
            }
            
            omp_unset_lock(&locks[(x+mx+M_x)%M_x]);
            ex*=e3x;
        }
    }
   
    //Spreading is the only part of the k-space sum that needs locks, so
    //get rid of them.
    for(int j = 0;j<M_x;j++)
        omp_destroy_lock(&locks[j]);
    delete locks;

    //---------------------------------------------------------------------
    //Step 2 : Frequency space filter
    //---------------------------------------------------------------------
    //We apply the appropriate k-space filter associated with the 
    //Stresslet. This involves Fast Fourier Transforms.
    
    //Call Matlab's routines to compute the 2D FFT. Other choices,
    //such as fftw could possibly be better, mainly because of the poor
    //complex data structure Matlab uses which we now have to deal with.
    mexCallMATLAB(1,&fft2lhs[0],1,&fft2rhs[0],"fft2");
    mexCallMATLAB(1,&fft2lhs[1],1,&fft2rhs[1],"fft2");
    mexCallMATLAB(1,&fft2lhs[2],1,&fft2rhs[2],"fft2");
    mexCallMATLAB(1,&fft2lhs[3],1,&fft2rhs[3],"fft2");

    //The output of the FFT is complex. Get pointers to the real and 
    //imaginary parts of the Hhats.
    double* Hhat1_re = mxGetPr(fft2lhs[0]);
    double* Hhat1_im = mxGetPi(fft2lhs[0]);
    
    double* Hhat2_re = mxGetPr(fft2lhs[1]);
    double* Hhat2_im = mxGetPi(fft2lhs[1]);
    
    double* Hhat3_re = mxGetPr(fft2lhs[2]);
    double* Hhat3_im = mxGetPi(fft2lhs[2]);
    
    double* Hhat4_re = mxGetPr(fft2lhs[3]);
    double* Hhat4_im = mxGetPi(fft2lhs[3]);

    //We cannot assume that the imaginary parts of the
    //Fourier transforms are non-zero. Let Matlab take care of the
    //memory management.
    int alloc3=0,alloc4=0;
    mwSize cs = M_x*M_y;
    
    if(Hhat1_im == NULL) {        
        Hhat1_im = (double*) mxCalloc(cs,sizeof(double));        
        mxSetPi(fft2lhs[0],Hhat1_im);
    }
    if(Hhat2_im == NULL) {
        Hhat2_im = (double*) mxCalloc(cs,sizeof(double));
        mxSetPi(fft2lhs[1],Hhat2_im);
    }
    
    if(Hhat3_im == NULL) {
        Hhat3_im = new double[M_x*M_y];
        memset(Hhat3_im,0,M_x*M_y*sizeof(double));
        alloc3 = 1;
    }
    
    if(Hhat4_im == NULL) {
        Hhat4_im = new double[M_x*M_y];
        memset(Hhat4_im,0,M_x*M_y*sizeof(double));
        alloc4 = 1;
    }

    //Apply filter in the frequency domain. This is a completely 
    //parallel operation. The FFT gives the frequency components in
    //non-sequential order, so we have to split the loops. One could
    //possibly use fast gaussian gridding here too, to get rid of the
    //exponentials, but as this loop accounts for a few percent of the
    //total runtime, this hardly seems worth the extra work.
    #pragma omp parallel for
    for(int j = 0;j<M_x;j++) {
        int ptr = j*M_y;
        double k1;
        if(j <= M_x/2)
            k1 = 2.0*pi/L_x*j;
        else
            k1 = 2.0*pi/L_x*(j-M_x);
        
        for(int k = 0;k<=M_y/2;k++,ptr++) {
            double k2 = 2.0*pi/L_y*k;
            double Ksq = k1*k1+k2*k2;
            
            double e = exp(-0.25*(1-eta)/(xi*xi)*Ksq)/Ksq/Ksq;
            double t = (1 + 0.25*Ksq/xi/xi);
           
            double B111 = 3.0*k1*Ksq - 2.0*k1*k1*k1;
            double B112 = k2*Ksq - 2.0*k1*k1*k2; //Also = B121
            double B122 = k1*Ksq - 2.0*k1*k2*k2;
            
            double tmp1re = Hhat1_re[ptr];
            double tmp1im = Hhat1_im[ptr];
            
            Hhat1_re[ptr] = e*t*(B111*Hhat1_re[ptr] + B112*(Hhat2_re[ptr]+Hhat3_re[ptr]) + B122*Hhat4_re[ptr]);
            Hhat1_im[ptr] = e*t*(B111*Hhat1_im[ptr] + B112*(Hhat2_im[ptr]+Hhat3_im[ptr]) + B122*Hhat4_im[ptr]);

            double B211 = k2*Ksq - 2.0*k2*k1*k1;
            double B212 = k1*Ksq - 2.0*k2*k1*k2; // Also = B221
            double B222 = 3.0*k2*Ksq - 2.0*k2*k2*k2;

            Hhat2_re[ptr] = e*t*(B211*tmp1re + B212*(Hhat2_re[ptr]+Hhat3_re[ptr]) + B222*Hhat4_re[ptr]);
            Hhat2_im[ptr] = e*t*(B211*tmp1im + B212*(Hhat2_im[ptr]+Hhat3_im[ptr]) + B222*Hhat4_im[ptr]);

        }
        for(int k = 0;k<M_y/2-1;k++,ptr++) {
            double k2 = 2.0*pi/L_y*(k-M_y/2+1);
            double Ksq = k1*k1+k2*k2;
           
            double e = exp(-0.25*(1-eta)/(xi*xi)*Ksq)/Ksq/Ksq;
            
            double t = (1 + 0.25*Ksq/xi/xi);
          
            double B111 = 3.0*k1*Ksq - 2.0*k1*k1*k1;
            double B112 = k2*Ksq - 2.0*k1*k1*k2; //Also = B121
            double B122 = k1*Ksq - 2.0*k1*k2*k2;
            
            double tmp1re = Hhat1_re[ptr];
            double tmp1im = Hhat1_im[ptr];
            
            Hhat1_re[ptr] = e*t*(B111*Hhat1_re[ptr] + B112*(Hhat2_re[ptr]+Hhat3_re[ptr]) + B122*Hhat4_re[ptr]);
            Hhat1_im[ptr] = e*t*(B111*Hhat1_im[ptr] + B112*(Hhat2_im[ptr]+Hhat3_im[ptr]) + B122*Hhat4_im[ptr]);

            double B211 = k2*Ksq - 2.0*k2*k1*k1;
            double B212 = k1*Ksq - 2.0*k2*k1*k2; // Also = B221
            double B222 = 3.0*k2*Ksq - 2.0*k2*k2*k2;

            Hhat2_re[ptr] = e*t*(B211*tmp1re + B212*(Hhat2_re[ptr]+Hhat3_re[ptr]) + B222*Hhat4_re[ptr]);
            Hhat2_im[ptr] = e*t*(B211*tmp1im + B212*(Hhat2_im[ptr]+Hhat3_im[ptr]) + B222*Hhat4_im[ptr]);

        }
    }
   
    //Remove the zero frequency term.
    Hhat1_re[0] = 0;
    Hhat2_re[0] = 0;
    Hhat1_im[0] = 0;
    Hhat2_im[0] = 0;
    
    if(alloc3 != 0) 
        delete Hhat3_im;
    if(alloc4 != 0) 
        delete Hhat4_im;
    
    //Get rid of the old H arrays. They are no longer needed.
    mxDestroyArray(fft2rhs[0]);
    mxDestroyArray(fft2rhs[1]);
    mxDestroyArray(fft2rhs[2]);
    mxDestroyArray(fft2rhs[3]);
    mxDestroyArray(fft2lhs[2]);
    mxDestroyArray(fft2lhs[3]);
    //Do the inverse 2D FFTs. We use Matlab's inbuilt functions again.
    mexCallMATLAB(1,&fft2rhs[0],1,&fft2lhs[0],"ifft2");
    mexCallMATLAB(1,&fft2rhs[1],1,&fft2lhs[1],"ifft2");
    
    //The pointer to the imaginary part. We don't need the real part.
    double* Ht1 = mxGetPi(fft2rhs[0]);
    double* Ht2 = mxGetPi(fft2rhs[1]);
    
    //In case the inverse FFTs are purely real. Not likely to happen, but
    //the program would crash without guarding for this.
    if(Ht1 == NULL) {
        Ht1 = new double[M_x*M_y];
        memset(Ht1,0,M_x*M_y*sizeof(double));
    }
    if(Ht2 == NULL) {
        Ht2 = new double[M_x*M_y];
        memset(Ht2,0,M_x*M_y*sizeof(double));
    }
    
    //Get rid of the Hhat arrays. They are no longer needed.
    mxDestroyArray(fft2lhs[0]);
    mxDestroyArray(fft2lhs[1]);
    
    //---------------------------------------------------------------------
    //Step 3 : Evaluating the velocity
    //---------------------------------------------------------------------
    //Here we compute the output velocity as a convolution of the Ht1
    //and Ht2 functions with a properly scaled gaussian. This is basically
    //gaussian blur, and we again use fast gaussian gridding. This 
    //procedure is fully parallel.
    
    //Create the output matrix.
    plhs[0] = mxCreateDoubleMatrix(2, Ntar, mxREAL);
    double* Tk = mxGetPr(plhs[0]);
    
    #pragma omp parallel for
    for(int k = 0;k<Ntar;k++) {
        
        //The Gaussian bells are translation invariant. We exploit this 
        //fact to avoid blow-up of the terms and the numerical instability
        //that follows. (px,py) is the center of the bell with the original
        //grid-alignment but close to the origin.
        double px = ptar[2*k] - h_x*floor(ptar[2*k]/h_x);
        double py = ptar[2*k+1] - h_y*floor(ptar[2*k+1]/h_y);
        
        double Lhalf_x = static_cast<double>(L_x/2.0);
        double Lhalf_y = static_cast<double>(L_y/2.0);
        
        //(mx,my) is the center grid point in H1 and H2.
        int mx;
        double tmpp = static_cast<double>((ptar[2*k]+Lhalf_x)/h_x);
        int tmppint = static_cast<int>(floor(tmpp));
        int tmppint2 = static_cast<int>(round(tmpp));
        int tmppint3 = static_cast<int>(ceil(tmpp-1));
        if(fabs(tmpp-tmppint2)<TOL) {
            tmppint = tmppint2;
            mx = static_cast<int>(tmppint-P/2);
        }else
        {
            if(fabs(px)>TOL) 
                mx = static_cast<int>(tmppint3-P/2);
            else
                mx = static_cast<int>(tmppint-P/2);           
        }
        
        int my;
        double tmppy = static_cast<double>((ptar[2*k+1]+Lhalf_y)/h_y);
        int tmppinty = static_cast<int>(floor(tmppy));
        int tmppint2y = static_cast<int>(round(tmppy));
        int tmppint3y = static_cast<int>(ceil(tmppy-1));
        if(fabs(tmppy-tmppint2y)<TOL) {
            tmppinty = tmppint2y;
            my = static_cast<int>(tmppinty-P/2);            
        }
        else
        {
            if(fabs(py)>TOL)
                my = static_cast<int>(tmppint3y-P/2);                
            else 
                my = static_cast<int>(tmppinty-P/2);
            
        }
        
        if (fabs(px - h_x) < TOL)
            px = 0;
        
        if (fabs(py - h_y) < TOL)
            py = 0;
        
        double tmp = -2*xi*xi/eta;
        double ex = exp(tmp*(px*px+py*py + 2*w*px));
        double e4y = exp(2*tmp*w*py);
        double e3x = exp(-2*tmp*h_x*px);
        double e3y = exp(-2*tmp*h_y*py);
        
        //If there is no wrap-around due to periodicity for this gaussian,
        //we use a faster loop. We go column by column, but as the 
        //matrices Ht1 and Ht2 are only read from we have no need for
        //locks and such.
        
        
        if(mx >= 0 && my >= 0 && mx < M_x-P-1 && my < M_y-P-1) {
            
            int idx = mx*M_y+my;
                
            for(int x = 0;x<P+1;x++) {
                double ey = ex*e4y*e1[x];
                
                for(int y = 0;y<P+1;y++) {
                    double tmp = ey*e1[y];
                    
                    Tk[2*k] += tmp*Ht1[idx];
                    Tk[2*k+1] += tmp*Ht2[idx++];
                    ey *= e3y;
                }
                ex*=e3x;
                idx += M_y-P-1;
            }
            
        }else{
            
            for(int x = 0;x<P+1;x++) {
                double ey = ex*e4y*e1[x];
                int xidx = ((x+mx+M_x)%M_x)*M_y;
                for(int y = 0;y<P+1;y++) {
                    double tmp = ey*e1[y];
                    int idx = ((y+my+M_y)%M_y)+xidx;
                    Tk[2*k] += tmp*Ht1[idx];
                    Tk[2*k+1] += tmp*Ht2[idx];
                    ey *= e3y;
                }
                ex*=e3x;
            }
        }
        
        tmp = 4*xi*xi/eta;
        tmp = tmp*tmp*h_x*h_y/pi;
        Tk[2*k] *= tmp / (4*pi);
        Tk[2*k+1] *= tmp / (4*pi);

    } 

    //Clean up
    mxDestroyArray(fft2rhs[0]);
    mxDestroyArray(fft2rhs[1]);
    delete e1;
}
