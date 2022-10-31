#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "base.h"
#include "atomic_data.h"



int new_zero_spec(double* spec){
    for (int i=0; i<NUM_ENERGY_BINS; i++)
        spec[i]=0;
    return 0;
}

int new_powerlaw_spec(double* spec, double Ndot_tot, double alpha_s){
    double nuHI=ion_potential(HI);
    double K=Ndot_tot*alpha_s*pow(nuHI,alpha_s);

    for (int i=0; i<NUM_ENERGY_BINS; i++)
        spec[i]=K*pow(eVarrayCen[i],-alpha_s-1);
    return 0;
}

double *init_spec() {
    double *spec = (double *)malloc(NUM_ENERGY_BINS*sizeof(double));
    assert(spec);
    new_zero_spec(spec);
    return spec;
}

int calc_transmitted_spec(double *Ndot_tot_out, double *Ndot_tot_in, double y[4], cellProp *cell){
    double tau_nu_HI, tau_nu_HeI, tau_nu_HeII, tau_nu_tot;
    for (int i=0; i<NUM_ENERGY_BINS; i++){
        tau_nu_HI   = y[0] * cell->nH * cell->dr_cm * sigma_HI_nu(eVarrayCen[i]);
        tau_nu_HeI  = y[1] * cell->nHe * cell->dr_cm * sigma_HeI_nu(eVarrayCen[i]);
        tau_nu_HeII = y[2] * cell->nHe * cell->dr_cm * sigma_HeII_nu(eVarrayCen[i]);
        tau_nu_tot = tau_nu_HI + tau_nu_HeI + tau_nu_HeII;
        Ndot_tot_out[i] = Ndot_tot_in[i] * exp(-tau_nu_tot);
    }
    return 0;
}
