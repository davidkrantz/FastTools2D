#include <mex.h>
#include <math.h>
#include <omp.h>
#include <pmmintrin.h>
#include <bessellib.h>
#include <gsl/gsl_sf_expint.h>
#include <ewald_tools.h>


#define _mm_shuf2_pd(__A) (static_cast<__m128d>(__builtin_ia32_shufpd (static_cast<__v2df>(__A), static_cast<__v2df>(__A), 1)))

inline double expint(double x);

/* Function to compute real space sum of K0(alpha*r) by the Spectral Ewald
 * method.
 *
 * Called in MATLAB as: G0r_ewald = mex_K0_real(psrc,ptar,fsrc,xi,nside,L,alpha),
 * where psrc and fsrc are 2xNsrc and ptar is 2xNtar.
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* List used for translating sources */
    static const int ilist_x[8] = {-1,-1,-1,0,0,1,1,1};
    static const int ilist_y[8] = {-1,0,1,-1,1,-1,0,1};
    
    /* Check inputs and outputs */
    if(nrhs != 9)
        mexErrMsgTxt("mex_K0_real: Incorrect number of input parameters");
    if(mxGetM(prhs[0]) != 2)
        mexErrMsgTxt("mex_K0_real: psrc must be a 2xn matrix.");
    if(mxGetM(prhs[1]) != 2)
        mexErrMsgTxt("mex_K0_real: ptar must be a 2xn matrix.");
    if(mxGetM(prhs[2]) != 1)
        mexErrMsgTxt("mex_K0_real: fsrc must be a 1xn matrix.");
    if(mxGetN(prhs[2]) != mxGetN(prhs[0]))
        mexErrMsgTxt("mex_K0_real: psrc and fsrc must be the same size.");
    
    /* (x,y)-coordinates of the source and target points psrc,ptar */
    double *psrc = mxGetPr(prhs[0]);
    double *ptar = mxGetPr(prhs[1]);
    /* Number of sources and targets */
    int Nsrc = mxGetN(prhs[0]);
    int Ntar = mxGetN(prhs[1]);
    /* (x,y)-components of the forces at points psrc */
    double *fsrc = mxGetPr(prhs[2]);
    /* The Ewald parameter xi */
    double xi = mxGetScalar(prhs[3]);
    /* Paramter s: #boxes along the periodic box */
    int nside_x = static_cast<int>(mxGetScalar(prhs[4]));
    int nside_y = static_cast<int>(mxGetScalar(prhs[5]));
    /* Length L of the periodic domain. NB assumed square */
    double Lx = mxGetScalar(prhs[6]);
    double Ly = mxGetScalar(prhs[7]);
    double alpha = mxGetScalar(prhs[8]);
    
    /* #Boxes contained in L^2, given by s*s */
    int num_boxes = nside_x*nside_y;
    
    /* Offsets refer to values that are added to a base pointer
     * in order to access a single element in a sequential list
     * of elements */
    int* particle_offsets_src = new int[Nsrc];
    int* box_offsets_src = new int[num_boxes+1]; //NB added +1 SP
    int* nsources_in_box = new int[num_boxes];
    
    int* particle_offsets_tar = new int[Ntar];
    int* box_offsets_tar = new int[num_boxes+1]; //NB added +1 SP
    int* ntargets_in_box = new int[num_boxes];
    
    /* Assigns particles to boxes on the current grid */
    Assign(psrc,ptar,Lx,Ly,Nsrc,Ntar,nside_x,nside_y,
            particle_offsets_src,box_offsets_src,nsources_in_box,
            particle_offsets_tar,box_offsets_tar,ntargets_in_box);
    
    
    /* 16-byte aligned arrays. Required to be compatible with SSE commands */
    double* psrc_a = static_cast<double*>(malloc (2*Nsrc*sizeof(double)));
    double* ptar_a = static_cast<double*>(malloc (2*Ntar*sizeof(double)));
    double* S  = static_cast<double*>(malloc (Nsrc*sizeof(double)));
    double* K0_R = static_cast<double*>(malloc (Ntar*sizeof(double)));
    
    /* Write from psrc to ps and from fsrc to fs, which are compatible with SSE commands */
    for(int j = 0;j<Nsrc;j++) {
        psrc_a[2*j] = psrc[2*particle_offsets_src[j]];
        psrc_a[2*j+1] = psrc[2*particle_offsets_src[j]+1];
        
        S[j] = fsrc[particle_offsets_src[j]]; 
    }
    
    for(int j = 0;j<Ntar;j++) {
        ptar_a[2*j] = ptar[2*particle_offsets_tar[j]];
        ptar_a[2*j+1] = ptar[2*particle_offsets_tar[j]+1];
        
        K0_R[j] = 0;
    }
    
    /* Cut off radius squared. nside such that all points in a
     * box are within the distance sqrt(cutoffsq) of each other */
    double cutoffsq = Lx*Ly/nside_x/nside_y;
    double xi2 = xi*xi;
    double a2 = alpha*alpha;
    
    double self = -0.5*E1(0.25*a2/xi2);
    
    /* Loop through boxes */
#pragma omp parallel for
    for(int current_box = 0;current_box<num_boxes;current_box++) {
        if(ntargets_in_box[current_box] == 0)
            continue;
        
        /* Temporary pointers to the particles of current box */
        int tidx = box_offsets_tar[current_box];
        int sidx = box_offsets_src[current_box];
        
        /* Compute the box self-interactions */
        for(int j=tidx;j<tidx+ntargets_in_box[current_box];j++) {
                        
            for(int k=sidx;k<sidx+nsources_in_box[current_box];k++) {
                                
                double x1 = ptar_a[2*j]-psrc_a[2*k];
                double x2 = ptar_a[2*j+1]-psrc_a[2*k+1];
                
                double r2 = x1*x1+x2*x2;
                                
                if(r2 == 0) {
                    K0_R[j] += self*S[k];
                    continue;
                } else {
                    /* Compute real space part GR = K0(xi^2 r^2, alpha^2/4xi^2) */
                    double x = xi2*r2;
                    double y = 0.25*a2/xi2;
                    int alg = 0;
                    int jmax = choose_algorithm(x, y, &alg);
                    double tmp = K0xy(x,y,alg,jmax);
//                     double tmp = 1;
                    K0_R[j] += 0.5*tmp*S[k];
                }
            }
        }
        
        /* Compute interactions from the nearest neighbors */
        for(int j=0;j<8;j++) {
            
            int per_source_x = current_box%nside_x+ilist_x[j];
            int per_source_y = current_box/nside_x+ilist_y[j];
            
            int t_x = (per_source_x+nside_x)%nside_x;
            int t_y = (per_source_y+nside_y)%nside_y;
            /* The number of the source nearest neighbor box */
            int source_box = t_y*nside_x + t_x;
                     
            if(nsources_in_box[source_box] > 0) {
                /* z-offset of the source box corrected for periodicity */
                double zoff_re = (Lx*(per_source_x-t_x))/nside_x;
                double zoff_im = (Ly*(per_source_y-t_y))/nside_y;
                
                for(int k=0;k<ntargets_in_box[current_box];k++) {
                    
                    int idx = box_offsets_src[source_box];
                    for(int l=0;l<nsources_in_box[source_box];l++,idx++) {
                        double x1 = ptar_a[2*(tidx+k)]-psrc_a[2*idx]-zoff_re;
                        double x2 = ptar_a[2*(tidx+k)+1]-psrc_a[2*idx+1]-zoff_im;
                        
                        double r2 = x1*x1+x2*x2;
                        if(r2 < cutoffsq) {
                            
                            /* Compute real space part GR = K0(xi^2 r^2, alpha^2/4xi^2) */
                            double x = xi2*r2;
                            double y = 0.25*a2/xi2;
                            int alg = 0;
                            int jmax = choose_algorithm(x, y, &alg);
                            double tmp = K0xy(x,y,alg,jmax);
                            K0_R[tidx+k] += 0.5*tmp*S[idx];                           
                        }                     
                    }                   
                }
            }        
        }  
    }
    
    
    free(psrc_a);
    free(ptar_a);
    free(S);
    
    plhs[0] = mxCreateDoubleMatrix(1, Ntar, mxREAL);
    
    double* T = mxGetPr(plhs[0]);
    
    for(int j = 0;j<Ntar;j++) {
        T[particle_offsets_tar[j]] = K0_R[j];
    }
    
    free(K0_R);
    
    delete particle_offsets_src;
    delete box_offsets_src;
    delete nsources_in_box;
    delete particle_offsets_tar;
    delete box_offsets_tar;
    delete ntargets_in_box;
    
}
