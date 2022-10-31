#ifndef _INCL_SPEC_
#define _INCL_SPEC_

#include<assert.h>
#include<math.h>
#include<stdlib.h>
#include"base.h"


int new_zero_spec(double* spec);

int new_powerlaw_spec(double* spec, double Ndot_tot, double alpha_s);

double *init_spec(); 

int calc_transmitted_spec(double *Ndot_tot_out, double *Ndot_tot_in, double y[4], cellProp *cell);
#endif
