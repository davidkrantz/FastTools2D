#include "mex.h"
#include <math.h>
#include <omp.h>
#include <pmmintrin.h>

#include "mm_mxmalloc.h"

#define pi 3.1415926535897932385

#define _mm_shuf2_pd(__A) (static_cast<__m128d>(__builtin_ia32_shufpd (static_cast<__v2df>(__A), static_cast<__v2df>(__A), 1)))

void Assign(double *psrc, double *ptar, double len_x, double len_y, int nsrc, 
        int ntar, int nside_x,int nside_y,
        int* particle_offsets_src,int* box_offsets_src,int* nsources_in_box,
        int* particle_offsets_tar,int* box_offsets_tar,int* ntargets_in_box);


inline double expint(double x);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /*List used for translating sources. FF */
    static const int ilist_x[8] = {-1,-1,-1,0,0,1,1,1};
    static const int ilist_y[8] = {-1,0,1,-1,1,-1,0,1};
    
    if(nrhs != 9)
        mexErrMsgTxt("Incorrect number of input parameters");
    
    if(mxGetM(prhs[0]) != 2)
        mexErrMsgTxt("psrc must be a 2xn matrix.");
    if(mxGetM(prhs[1]) != 2)
        mexErrMsgTxt("ptar must be a 2xn matrix.");
    if(mxGetM(prhs[2]) != 2)
        mexErrMsgTxt("f must be a 2xn matrix.");
    if(mxGetM(prhs[3]) != 2)
        mexErrMsgTxt("n must be a 2xn matrix.");
    if(mxGetN(prhs[2]) != mxGetN(prhs[0]))
        mexErrMsgTxt("psrc and f must be the same size.");
    if(mxGetN(prhs[3]) != mxGetN(prhs[0]))
        mexErrMsgTxt("psrc and n must be the same size.");
    
    /*(x,y)-coordinates of the source and target points psrc,ptar. SP*/
    double *psrc = mxGetPr(prhs[0]);
    double *ptar = mxGetPr(prhs[1]);
    /*Number of sources and targets. FF*/
    int Nsrc = mxGetN(prhs[0]);
    int Ntar = mxGetN(prhs[1]);
    /*(x,y)-components of the forces at points p. FF*/
    double *f = mxGetPr(prhs[2]);
    double *n = mxGetPr(prhs[3]);
    /*The Ewald parameter xi. FF*/
    double xi = mxGetScalar(prhs[4]);
    /*Paramter s: #boxes along the periodic box. FF*/
    int nside_x = static_cast<int>(mxGetScalar(prhs[5]));
    int nside_y = static_cast<int>(mxGetScalar(prhs[6]));
    /*Length L of the periodic domain. FF*/
    double len_x = mxGetScalar(prhs[7]);
    double len_y = mxGetScalar(prhs[8]);
    
    /*#Boxes contained in L^2, given by s*s. FF*/
    int num_boxes = nside_x*nside_y;
    
    /*offsets refer to values that are added to a base pointer
     * in order to access a single element in a sequential list
     * of elements. FF */
    int* particle_offsets_src = new int[Nsrc];
    int* box_offsets_src = new int[num_boxes+1]; //NB added +1 SP
    int* nsources_in_box = new int[num_boxes];
    
    int* particle_offsets_tar = new int[Ntar];
    int* box_offsets_tar = new int[num_boxes+1]; //NB added +1 SP
    int* ntargets_in_box = new int[num_boxes];
    
    /*Assigns particles to boxes on the current grid. FF*/
    Assign(psrc,ptar,len_x,len_y,Nsrc,Ntar,nside_x,nside_y,
            particle_offsets_src,box_offsets_src,nsources_in_box,
            particle_offsets_tar,box_offsets_tar,ntargets_in_box);
    
    
    /*16-byte aligned arrays. Required to be compatible with SSE commands. FF*/
    double* psrc_a = static_cast<double*>(_mm_mxMalloc (2*Nsrc*sizeof(double), 16));
    double* ptar_a = static_cast<double*>(_mm_mxMalloc (2*Ntar*sizeof(double), 16));
    double* S  = static_cast<double*>(_mm_mxMalloc (4*Nsrc*sizeof(double), 16));
    double* Ts = static_cast<double*>(_mm_mxCalloc (2*Ntar, sizeof(double),16));
    
    /*Write from p to ps and from f to fs, which are compatible with SSE commands. FF*/
    for(int j = 0;j<Nsrc;j++) {
        
        psrc_a[2*j] = psrc[2*particle_offsets_src[j]];
        psrc_a[2*j+1] = psrc[2*particle_offsets_src[j]+1];
        
        S[4*j] = f[2*particle_offsets_src[j]]*n[2*particle_offsets_src[j]];
        S[4*j+1] = f[2*particle_offsets_src[j]]*n[2*particle_offsets_src[j]+1];
        S[4*j+2] = f[2*particle_offsets_src[j]+1]*n[2*particle_offsets_src[j]];
        S[4*j+3] = f[2*particle_offsets_src[j]+1]*n[2*particle_offsets_src[j]+1];
        
    }
    
    for(int j = 0;j<Ntar;j++) {
        ptar_a[2*j] = ptar[2*particle_offsets_tar[j]];
        ptar_a[2*j+1] = ptar[2*particle_offsets_tar[j]+1];
    }
    /*Cut off radius squared. nside such that all points in a
     * box are within the distance sqrt(cutoffsq) of each other. FF*/
    double cutoffsq = len_x*len_y/nside_x/nside_y;
    double xi2 = xi*xi;
    
   // mexPrintf("rc2 = %3.3f, xi = %3.3f\n", cutoffsq, xi);

    /*Loop through boxes*/
#pragma omp parallel for    
    for(int current_box = 0;current_box<num_boxes;current_box++) {
        if(ntargets_in_box[current_box] == 0)
            continue;
        
        //Temporary pointers to the particles of current box.
        int tidx = box_offsets_tar[current_box];
        int sidx = box_offsets_src[current_box];
                
        //Compute the box self-interactions.
        for(int j=tidx;j<tidx+ntargets_in_box[current_box];j++) {
            
            for(int k=sidx;k<sidx+nsources_in_box[current_box];k++) {
                                
//                 double x1 = (psrc_a[2*k]-ptar_a[2*j]);
//                 double x2 = (psrc_a[2*k+1]-ptar_a[2*j+1]);
                double x1 = ptar_a[2*j]-psrc_a[2*k];
                double x2 = ptar_a[2*j+1]-psrc_a[2*k+1];
                
                double r2 = x1*x1+x2*x2;
                
                if(r2 == 0)
                    continue;
          
                double e2 = exp(-xi2*r2);
                double prefac = 2*xi2;
                double facb = -4*(1+xi2*r2)/r2/r2;
                                
                double T111 = x1*x1*x1*facb + prefac*3*x1;
                double T112 = x1*x1*x2*facb + prefac*x2;
                double T122 = x1*x2*x2*facb + prefac*x1;
                
                Ts[2*j] += e2*(T111*S[4*k] + T112*(S[4*k+1] + S[4*k+2]) + T122*S[4*k+3]);
                                
                double T211 = x2*x1*x1*facb + prefac*x2;
                double T212 = x2*x1*x2*facb + prefac*x1;
                double T222 = x2*x2*x2*facb + prefac*3*x2;
                Ts[2*j+1] += e2*(T211*S[4*k] + T212*(S[4*k+1] + S[4*k+2]) + T222*S[4*k+3]); 
                
               // printf("x1 = %3.5e, x2 = %3.5e, r=%3.5e\n", x1, x2, r2); 
           
            } 
        }
        
        //Compute interactions from the nearest neighbors.
        for(int j=0;j<8;j++) {
            
            int per_source_x = current_box%nside_x+ilist_x[j];
            int per_source_y = current_box/nside_x+ilist_y[j];
            
            int t_x = (per_source_x+nside_x)%nside_x;
            int t_y = (per_source_y+nside_y)%nside_y;
            //The number of the source nearest neighbor box.
            int source_box = t_y*nside_x + t_x;
    
            if(nsources_in_box[source_box] > 0) {
                //z-offset of the source box corrected for periodicity.
                double zoff_re = (len_x*(per_source_x-t_x))/nside_x;
                double zoff_im = (len_y*(per_source_y-t_y))/nside_y;
                
                for(int k=0;k<ntargets_in_box[current_box];k++) {
                    
                    int idx = box_offsets_src[source_box];
                    for(int l=0;l<nsources_in_box[source_box];l++,idx++) {
//                         double x1 = psrc_a[2*idx]-ptar_a[2*(tidx+k)]+zoff_re;
//                         double x2 = psrc_a[2*idx+1]-ptar_a[2*(tidx+k)+1]+zoff_im;
                        double x1 = ptar_a[2*(tidx+k)]-psrc_a[2*idx]-zoff_re;
                        double x2 = ptar_a[2*(tidx+k)+1]-psrc_a[2*idx+1]-zoff_im;
                        
                        double r2 = x1*x1+x2*x2;
                        //printf("x1 = %3.5e, x2 = %3.5e, r=%3.5e\n", x1, x2, r2);
                        
                        if(r2 < cutoffsq) {

                            double e2 = exp(-xi2*r2);
                            double prefac = 2*xi2;
                            double facb = -4*(1+xi2*r2)/r2/r2;
                            
                            double T111 = x1*x1*x1*facb + prefac*3*x1;
                            double T112 = x1*x1*x2*facb + prefac*x2;
                            double T122 = x1*x2*x2*facb + prefac*x1;
                            Ts[2*(tidx+k)] += e2*(T111*S[4*idx] + T112*(S[4*idx+1] + S[4*idx+2]) + T122*S[4*idx+3]);
                            
                            double T211 = x2*x1*x1*facb + prefac*x2;
                            double T212 = x2*x1*x2*facb + prefac*x1;
                            double T222 = x2*x2*x2*facb + prefac*3*x2;
                            Ts[2*(tidx+k)+1] += e2*(T211*S[4*idx] + T212*(S[4*idx+1] + S[4*idx+2]) + T222*S[4*idx+3]);
                        }                        
                    }                  
                }
            } 
        }
    }
    
    _mm_mxFree(psrc_a);
    _mm_mxFree(ptar_a);
    _mm_mxFree(S);
    
    plhs[0] = mxCreateDoubleMatrix(2, Ntar, mxREAL);
    
    double* T = mxGetPr(plhs[0]);
    double A = len_x*len_y;
    
    double s1 = 0;
    double s2 = 0;
    for(int j = 0;j<Nsrc;j++) {
        double t = (n[2*j]*f[2*j]+n[2*j+1]*f[2*j+1]);
        s1 += t*psrc[2*j];
        s2 += t*psrc[2*j+1];
        
       // mexPrintf("x = %3.3f, y=%3.3f\n", psrc[2*j], psrc[3*j+1]);
        
    }
    
    s1/=A;
    s2/=A;
    
   // mexPrintf("s1 = %3.3g, s2 = %3.3g, A = %3.3f\n", s1, s2, A);
    
    for(int j = 0;j<Ntar;j++) {
        T[2*particle_offsets_tar[j]] = Ts[2*j]/4/pi+s1;
        T[2*particle_offsets_tar[j]+1] = Ts[2*j+1]/4/pi+s2;
    }
    
    _mm_mxFree(Ts);
    
    delete particle_offsets_src;
    delete box_offsets_src;
    delete nsources_in_box;
    delete particle_offsets_tar;
    delete box_offsets_tar;
    delete ntargets_in_box;
    
}

/*------------------------------------------------------------------------
 *This function assigns particles to boxes on the current grid.
 *------------------------------------------------------------------------
 */
void Assign(double *psrc, double *ptar, double len_x, double len_y, int nsrc, 
        int ntar, int nside_x, int nside_y,
        int* particle_offsets_src,int* box_offsets_src,int* nsources_in_box,
        int* particle_offsets_tar,int* box_offsets_tar,int* ntargets_in_box){
    
    //The total number of boxes.
    int number_of_boxes = nside_x*nside_y;
    int j;
    int* in_box_src = new int[nsrc];
    int* in_box_tar = new int[ntar];
    int* particlenum_in_box_src = new int[number_of_boxes+1];
    int* particlenum_in_box_tar = new int[number_of_boxes+1];
    
    //Clear the nparticles_in_box arrays.
    for(j = 0;j<number_of_boxes;j++){
        nsources_in_box[j] = 0;
        ntargets_in_box[j] = 0;
    }
    
    //Assign the sources to boxes.
    for(j = 0;j<nsrc;j++) {
        int box_x = (int)floor(nside_x*(psrc[2*j]/len_x+0.5));
        int box_y = (int)floor(nside_y*(psrc[2*j+1]/len_y+0.5));
        if(box_x < 0) box_x = 0;
        if(box_x >= nside_x) box_x = nside_x-1;
        if(box_y < 0) box_y = 0;
        if(box_y >= nside_y) box_y = nside_y-1;
        in_box_src[j] =  box_y*nside_x + box_x;
        nsources_in_box[in_box_src[j]]++;
    }
    
    //Assign the targets to boxes.
    for(j = 0;j<ntar;j++) {
        int box_x = (int)floor(nside_x*(ptar[2*j]/len_x+0.5));
        int box_y = (int)floor(nside_y*(ptar[2*j+1]/len_y+0.5));
        if(box_x < 0) box_x = 0;
        if(box_x >= nside_x) box_x = nside_x-1;
        if(box_y < 0) box_y = 0;
        if(box_y >= nside_y) box_y = nside_y-1;
        in_box_tar[j] =  box_y*nside_x + box_x;
        ntargets_in_box[in_box_tar[j]]++;
    }
    
    /*
    printf("SOURCES\n");
    for (j=0; j<nsrc; j++)
    {
     printf("%d\n", in_box_src[j]);  
    }
    
    printf("\nTargets\n");
    for (j=0; j<ntar; j++)
    {
     printf("%d\n", in_box_tar[j]);  
    }
    
    printf("\nBOX COUNT\n");
    for (j=0; j<number_of_boxes; j++)
    {
     printf("%d\n", ntargets_in_box[j]);  
     printf("%d\n", nsources_in_box[j]);  
    }
    */
     
    //box_offsets is the offsets of the boxes in particle_offsets to be
    //computed below. particlenum_in_box is a temporary array.
    box_offsets_src[0] = particlenum_in_box_src[0] = 0;
    box_offsets_tar[0] = particlenum_in_box_tar[0] = 0;
    int isrc = 0;
    int itar = 0;
    
    for(j=1;j<number_of_boxes+1;j++) {
        isrc += nsources_in_box[j-1];
        
        particlenum_in_box_src[j] = isrc;
        
        box_offsets_src[j] = isrc;
        
        itar += ntargets_in_box[j-1];
        particlenum_in_box_tar[j] = itar;
        
        box_offsets_tar[j] = itar;
    }
    
    
    //The offsets of the particles in each box in the z and q arrays.
    for(j=0;j<nsrc;j++)
        particle_offsets_src[particlenum_in_box_src[in_box_src[j]]++] = j;
    for(j=0;j<ntar;j++)
        particle_offsets_tar[particlenum_in_box_tar[in_box_tar[j]]++] = j;
    
    delete in_box_src;
    delete particlenum_in_box_src;
    delete in_box_tar;
    delete particlenum_in_box_tar;
    
}
