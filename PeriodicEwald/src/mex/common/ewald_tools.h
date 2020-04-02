#ifndef EWALD_TOOLS  
#define EWALD_TOOLS

#include <math.h>

void Assign(double *psrc, double *ptar, double len_x, double len_y, int nsrc, 
        int ntar, int nside_x, int nside_y, int* particle_offsets_src,
        int* box_offsets_src,int* nsources_in_box, int* particle_offsets_tar,
        int* box_offsets_tar,int* ntargets_in_box);

#endif
