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
