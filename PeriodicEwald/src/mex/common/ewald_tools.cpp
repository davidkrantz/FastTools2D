#include "ewald_tools.h"

/*------------------------------------------------------------------------
 *This function assigns particles to boxes on the current grid.
 *------------------------------------------------------------------------
 */
void Assign(double *psrc, double *ptar, double Lx, double Ly, int nsrc, 
        int ntar, int nside_x, int nside_y, int* particle_offsets_src,
        int* box_offsets_src,int* nsources_in_box, int* particle_offsets_tar,
        int* box_offsets_tar,int* ntargets_in_box){
    
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
        int box_x = (int)floor(nside_x*(psrc[2*j]/Lx+0.5));
        int box_y = (int)floor(nside_y*(psrc[2*j+1]/Ly+0.5));
        if(box_x < 0) box_x = 0;
        if(box_x >= nside_x) box_x = nside_x-1;
        if(box_y < 0) box_y = 0;
        if(box_y >= nside_y) box_y = nside_y-1;
        in_box_src[j] =  box_y*nside_x + box_x;
        nsources_in_box[in_box_src[j]]++;
    }
    
    //Assign the targets to boxes.
    for(j = 0;j<ntar;j++) {
        int box_x = (int)floor(nside_x*(ptar[2*j]/Lx+0.5));
        int box_y = (int)floor(nside_y*(ptar[2*j+1]/Ly+0.5));
        if(box_x < 0) box_x = 0;
        if(box_x >= nside_x) box_x = nside_x-1;
        if(box_y < 0) box_y = 0;
        if(box_y >= nside_y) box_y = nside_y-1;
        in_box_tar[j] =  box_y*nside_x + box_x;
        ntargets_in_box[in_box_tar[j]]++;
    }
    
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

/*------------------------------------------------------------------------
 *This function finds the node to begin the Gaussian blur
 *------------------------------------------------------------------------
 */
void FindClosestNode(double x, double y, double Lx, double Ly, double h, int P, 
			int* mx, int* my, double* px, double* py){
	
	double TOL = 1e-13;

        *px = x - h*floor(x/h);
        *py = y - h*floor(y/h);

        double Lhalf_x = static_cast<double>(Lx/2.0);
        double Lhalf_y = static_cast<double>(Ly/2.0);
        
        //(mx,my) is the center grid point in H1 and H2.
        *mx = static_cast<int>(floor((x+Lhalf_x)/h) - P/2);
        *my = static_cast<int>(floor((y+Lhalf_y)/h) - P/2);

	// correct for cases where target is very close to a grid node

	if (fabs(*px) < TOL || fabs(*px - h) < TOL)
            *px = 0;
        
        if (fabs(*py) < TOL || fabs(*py - h) < TOL)
            *py = 0;

	if (*px == 0 && remainder((x+Lhalf_x)/h - P/2, 1) < 0)
		(*mx)++;

	if (*py == 0 && remainder((y+Lhalf_y)/h - P/2, 1) < 0)
		(*my)++;
}
