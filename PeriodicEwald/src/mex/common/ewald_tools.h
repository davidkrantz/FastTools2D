#ifndef EWALD_TOOLS  
#define EWALD_TOOLS

#include <math.h>
#include <string.h>
#include "mex.h"

void Assign(double *psrc, double *ptar, double len_x, double len_y, int nsrc, 
        int ntar, int nside_x, int nside_y, int* particle_offsets_src,
        int* box_offsets_src,int* nsources_in_box, int* particle_offsets_tar,
        int* box_offsets_tar,int* ntargets_in_box);

void FindClosestNode(double x, double y, double Lx, double Ly, double h, int P, 
			int* mx, int* my, double* px, double* py);
#endif
