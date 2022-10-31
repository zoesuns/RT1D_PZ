#ifndef _INCL_FAKE_CELLS_
#define _INCL_FAKE_CELLS_

#include <stdlib.h>
#include <math.h>
#include "../base.h"
void H_only_cell(cellProp *cell, double dist_pMpc, double dr_pMpc);
void He_only_cell(cellProp *cell, double dist_pMpc, double dr_pMpc);
void cold_dense_neutral_cell(cellProp *cell, double dist_pMpc, double dr_pMpc);
void hot_diffuse_ionized_cell(cellProp *cell, double dist_pMpc, double dr_pMpc);
void z7_cell(cellProp *cell, double dist_pMpc, double dr_pMpc);

#endif
