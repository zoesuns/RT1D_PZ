#include<stdlib.h>
#include<math.h>
#include"base.h"

void fake_cell(cellProp *cell, double dist_pMpc, double dr_pMpc){
    cell->dist_pMpc=dist_pMpc;
    cell->dr_pMpc=dr_pMpc;
    cell->nH=0.1;
    //cell->nHe=cell->nH*0.0789474;
    cell->nHe=1e-8;
    cell->xHI=1;
    cell->xHeI=1;
    cell->xHeII=0;
    cell->T=1e3;
}
