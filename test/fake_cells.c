#include <stdlib.h>
#include <math.h>
#include "../base.h"
#include "../cosmology.h"

void z7_cell(cellProp *cell, double dist_pMpc, double dr_pMpc){
    double rhob=cosmo_rhoc_z(7)*cosmo_Ob_z(7);
    cell->dist_pMpc=dist_pMpc;
    cell->dr_pMpc=dr_pMpc;
    cell->nH=rhob*(1.-cosmo_Y)/protonMass;
    cell->nHe=rhob*cosmo_Y/protonMass/4.;
    cell->xHI=1;
    cell->xHeI=1;
    cell->xHeII=0;
    cell->T=10;
    cell_parse_info(cell);
}

void H_only_cell(cellProp *cell, double dist_pMpc, double dr_pMpc){
    cell->dist_pMpc=dist_pMpc;
    cell->dr_pMpc=dr_pMpc;
    cell->nH=1e-4;
    cell->nHe=1e-8;
    cell->xHI=1;
    cell->xHeI=1;
    cell->xHeII=0;
    cell->T=1e2;
    cell_parse_info(cell);
}

void He_only_cell(cellProp *cell, double dist_pMpc, double dr_pMpc){
    cell->dist_pMpc=dist_pMpc;
    cell->dr_pMpc=dr_pMpc;
    cell->nH=1e-8;
    cell->nHe=1.e-3;
    cell->xHI=1;
    cell->xHeI=1;
    cell->xHeII=0;
    cell->T=1e2;
    cell_parse_info(cell);
}


void cold_dense_neutral_cell(cellProp *cell, double dist_pMpc, double dr_pMpc){
    cell->dist_pMpc=dist_pMpc;
    cell->dr_pMpc=dr_pMpc;
    cell->nH=1;
    cell->nHe=cell->nH*0.0789474;
    cell->xHI=1;
    cell->xHeI=1;
    cell->xHeII=0;
    cell->T=1e2;
    cell_parse_info(cell);
}

void hot_diffuse_ionized_cell(cellProp *cell, double dist_pMpc, double dr_pMpc){
    cell->dist_pMpc=dist_pMpc;
    cell->dr_pMpc=dr_pMpc;
    cell->nH=1e-4;
    cell->nHe=cell->nH*0.0789474;
    cell->xHI=1e-3;
    cell->xHeI=1e-3;
    cell->xHeII=1e-1;
    cell->T=1e5;
    cell_parse_info(cell);
}
