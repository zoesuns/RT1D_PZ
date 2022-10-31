#ifndef _INCL_ATOMIC_DATA_
#define _INCL_ATOMIC_DATA_

#include <math.h>
#include "base.h"

#define threshHIeV   (13.6)
#define threshHeIeV  (24.6)
#define threshHeIIeV (54.4)

typedef enum _ionName {HI,HeI,HeII,HeI_r,HeI_d} ionName;

double ion_potential(ionName species); // should be deleted later

double sigma_HI_nu(double nu);
double sigma_HeI_nu(double nu);
double sigma_HeII_nu(double nu);

double eGamma_HI(double T);
double eGamma_HeI(double T);
double eGamma_HeII(double T);

double brem_cooling(double T, double ne, double nHII, double nHeII, double nHeIII);

double inv_comp_cooling(double T, double ne, double zuni);


int init_load_recomb_data();
int free_loaded_recomb_data();

double recomb_rate(double T, ionName species);

double recomb_cooling(double T, double ne, double nHII, double nHeII, double nHeIII);

double colli_ioniz_cooling(double T, double ne, double nHI, double nHeI, double nHeII);

double colli_excit_cooling(double T, double ne, double nHI, double nHeI, double nHeII);

double brem_cooling(double T, double ne, double nHII, double nHeII, double nHeIII);

double inv_comp_cooling(double T, double ne, double zuni);


double calc_Cooling(cellProp *cell, double zuni);
#endif
