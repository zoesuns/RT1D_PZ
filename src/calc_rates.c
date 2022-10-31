#include <stdio.h>
#include <math.h>
#include "atomic_data.h"
#include "base.h"

inline void scale_vector(double *out, size_t n, double a, double *b){
    for (int i=0; i<n; i++)
        out[i]=a*b[i];
}

inline void add_vector(double *out, size_t n, double *a, double *b){
    for (int i=0; i<n; i++)
        out[i]=a[i]+b[i];
}

inline void sub_vector(double *out, size_t n, double *a, double *b){
    for (int i=0; i<n; i++)
        out[i]=a[i]-b[i];
}

void exp_vector(double *out, size_t n, double *b){
    for (int i=0; i<n; i++)
        out[i]=exp(b[i]);
}



inline double inner_product(size_t n, double *a, double *b){
    double sum=0;
    for (int i=0; i<n; i++)
        sum+=(a[i]*b[i]);
    return sum;
}

int calc_photoGH(double absNdot_nu[], cellProp *cell, double Gamma[], double phHeating[]){
    
    double tau_nu_HI, tau_nu_HeI, tau_nu_HeII, tau_nu_tot;
    double q_nu_HI, q_nu_HeI, q_nu_HeII, q_nu_tot;
    double p_nu_HI, p_nu_HeI, p_nu_HeII, p_nu_tot;

    double pHIqq_nu, pHeIqq_nu, pHeIIqq_nu, D_nu;
    double Pabs_nu_HI, Pabs_nu_HeI, Pabs_nu_HeII;

    double NPdnu_HI, NPdnu_HeI, NPdnu_HeII;

    double int_NPdnu_HI=0, int_NPdnu_HeI=0, int_NPdnu_HeII=0;
    double int_EnNPdnu_HI=0, int_EnNPdnu_HeI=0, int_EnNPdnu_HeII=0;

    int i=0;
    //double testN=0;

    for (i=0; i<NUM_ENERGY_BINS; i++){
        tau_nu_HI   = cell->nHI   * cell->dr_cm * sigma_HI_nu(eVarrayCen[i]);
        tau_nu_HeI  = cell->nHeI  * cell->dr_cm * sigma_HeI_nu(eVarrayCen[i]);
        tau_nu_HeII = cell->nHeII * cell->dr_cm * sigma_HeII_nu(eVarrayCen[i]);
        tau_nu_tot = tau_nu_HI + tau_nu_HeI + tau_nu_HeII;
        //printf("tau_nu_HI=%e, tau_nu_HeI=%e, tau_nu_HeII=%e, tau_nu_tot=%e\n", tau_nu_HI, tau_nu_HeI, tau_nu_HeII, tau_nu_tot);

        //printf("tau_nu_tot=%e\n", tau_nu_tot);
        if (tau_nu_tot>1e-9){
            //printf("yes\n");	
	        q_nu_HI=exp(-tau_nu_HI);
	        q_nu_HeI=exp(-tau_nu_HeI);
	        q_nu_HeII=exp(-tau_nu_HeII);
	        q_nu_tot=exp(-tau_nu_tot);
	
	        p_nu_HI   = 1-q_nu_HI;
	        p_nu_HeI  = 1-q_nu_HeI;
	        p_nu_HeII = 1-q_nu_HeII;
	        p_nu_tot  = 1-q_nu_tot;
            //printf("q_nu_HI=%e, q_nu_HeI=%e, q_nu_HeII=%e\n", q_nu_HI, q_nu_HeI, q_nu_HeII);
            //printf("p_nu_HI=%e, p_nu_HeI=%e, p_nu_HeII=%e\n", p_nu_HI, p_nu_HeI, p_nu_HeII);
            //printf("q_nu_HI=%e, p_nu_HI=%e\n", q_nu_HI, p_nu_HI);
            //printf("q_nu_HeI=%e, p_nu_HeI=%e\n", q_nu_HeI, p_nu_HeI);
            //printf("q_nu_HeII=%e, p_nu_HeII=%e\n", q_nu_HeII, p_nu_HeII);
            //printf("q_nu_tot=%e\n", q_nu_tot);
            //printf("pqq/D=%e\n", pHIqq_nu/D_nu);

	        pHIqq_nu   = p_nu_HI*q_nu_HeI*q_nu_HeII;
	        pHeIqq_nu  = q_nu_HI*p_nu_HeI*q_nu_HeII;
	        pHeIIqq_nu = q_nu_HI*q_nu_HeI*p_nu_HeII;
	        D_nu = pHIqq_nu + pHeIqq_nu + pHeIIqq_nu;
            //printf("D_nu=%e\n", D_nu);
            //printf("D_nu=%e\n", D_nu);
            //printf("pHIqq_nu=%e, pHeIqq_nu=%e, pHeIIqq_nu=%e\n", pHIqq_nu, pHeIqq_nu, pHeIIqq_nu);
            if (D_nu>0){
	            Pabs_nu_HI   =   pHIqq_nu*p_nu_tot/D_nu;
	            Pabs_nu_HeI  =  pHeIqq_nu*p_nu_tot/D_nu;
	            Pabs_nu_HeII = pHeIIqq_nu*p_nu_tot/D_nu;
            } else {
                if (q_nu_HI<1e-9 && q_nu_HeI<1e-9 && q_nu_HeII<1e-9){
                    Pabs_nu_HI   = (tau_nu_HI>tau_nu_HeI && tau_nu_HI>tau_nu_HeII) ? 1 : 0;
                    Pabs_nu_HeI  = (tau_nu_HeI>tau_nu_HI && tau_nu_HeI>tau_nu_HeII) ? 1 : 0;
                    Pabs_nu_HeII = 1 - Pabs_nu_HI - Pabs_nu_HeI;
                }
                if (q_nu_HI<1e-9 && q_nu_HeI<1e-9){
                    Pabs_nu_HI   = tau_nu_HI>tau_nu_HeI? 1:0;
                    Pabs_nu_HeI  = 1-Pabs_nu_HI;
                    Pabs_nu_HeII = 0;
                } else if (q_nu_HI<1e-9 && q_nu_HeII<1e-9){
                    Pabs_nu_HI   = tau_nu_HI>tau_nu_HeII? 1 : 0;
                    Pabs_nu_HeI  = 0;
                    Pabs_nu_HeII = 1 - Pabs_nu_HI;
                } else if (q_nu_HeI<1e-9 && q_nu_HeII<1e-9){
                    Pabs_nu_HI   = 0;
                    Pabs_nu_HeI  = tau_nu_HeI>tau_nu_HeII? 1 : 0;
                    Pabs_nu_HeII = 1 - Pabs_nu_HeI;
                } else {
            //printf("q_nu_HI=%e, q_nu_HeI=%e, q_nu_HeII=%e\n", q_nu_HI, q_nu_HeI, q_nu_HeII);
            //printf("p_nu_HI=%e, p_nu_HeI=%e, p_nu_HeII=%e, D_nu=%e\n", p_nu_HI, p_nu_HeI, p_nu_HeII, D_nu);
            //printf("pHIqq_nu=%e, pHeIqq_nu=%e, pHeIIqq_nu=%e\n", pHIqq_nu, pHeIqq_nu, pHeIIqq_nu);
                    printf("ERROR: Pabs shouldn't be this value!!!\n");
                    exit(VALUE_ERROR);
                }


            }            
            //printf("Pabs_nu_HI=%e\n", Pabs_nu_HI);
            //printf("opt thin=%e\n", tau_nu_HI);
        }
        else{
            //printf("no\n");
            Pabs_nu_HI=tau_nu_HI;
            Pabs_nu_HeI=tau_nu_HeI;
            Pabs_nu_HeII=tau_nu_HeII;
        }
        //testN += absNdot_nu[i]*eVarrayWidth[i];
        //printf("Pabs_nu_HI = %e, Pabs_nu_HeI = %e, Pabs_nu_HeII = %e\n", Pabs_nu_HI, Pabs_nu_HeI, Pabs_nu_HeII);
        NPdnu_HI = absNdot_nu[i]*Pabs_nu_HI  *eVarrayWidth[i];
        //printf("NPdnu_HI = %e\n", NPdnu_HI);
        NPdnu_HeI = absNdot_nu[i]*Pabs_nu_HeI *eVarrayWidth[i];
        NPdnu_HeII = absNdot_nu[i]*Pabs_nu_HeII*eVarrayWidth[i];
        int_NPdnu_HI  += NPdnu_HI;
        int_NPdnu_HeI += NPdnu_HeI;
        int_NPdnu_HeII+= NPdnu_HeII;
        int_EnNPdnu_HI   += (NPdnu_HI  * (eVarrayCen[i]-threshHIeV  )*eV2erg);
        int_EnNPdnu_HeI  += (NPdnu_HeI * (eVarrayCen[i]-threshHeIeV )*eV2erg);
        int_EnNPdnu_HeII += (NPdnu_HeII* (eVarrayCen[i]-threshHeIIeV)*eV2erg);
    }
    //printf("testN=%10.4e\n",testN); 
    if (cell->nHI<=0){
        Gamma[0] = 0;
        phHeating[0] = 0;
    } else {
        Gamma[0] =   int_NPdnu_HI /(4*M_PI*cell->dist_cm*cell->dist_cm*cell->dr_cm)/cell->nHI;
        phHeating[0] =   int_EnNPdnu_HI /(4*M_PI*cell->dist_cm*cell->dist_cm*cell->dr_cm)/cell->nHI;
    }

    if (cell->nHeI<=0){
        Gamma[1] = 0;
        phHeating[1] = 0;
    } else {
        Gamma[1] =  int_NPdnu_HeI /(4*M_PI*cell->dist_cm*cell->dist_cm*cell->dr_cm)/cell->nHeI;
        phHeating[1] =  int_EnNPdnu_HeI /(4*M_PI*cell->dist_cm*cell->dist_cm*cell->dr_cm)/cell->nHeI;
    }

    if (cell->nHeII<=0){
        Gamma[2] = 0;
        phHeating[2] = 0;
    } else {
        Gamma[2] = int_NPdnu_HeII /(4*M_PI*cell->dist_cm*cell->dist_cm*cell->dr_cm)/cell->nHeII;
        phHeating[2] = int_EnNPdnu_HeII /(4*M_PI*cell->dist_cm*cell->dist_cm*cell->dr_cm)/cell->nHeII;
    }
    return 0;
}


int est_ion_time(double *time, double *absNdot_nu, double dist_cm){
    double tmp0=0, tmp1=0, tmp2=0;
    for (int i=0; i<NUM_ENERGY_BINS; i++){
        tmp0+=absNdot_nu[i]*sigma_HI_nu(eVarrayCen[i])*eVarrayWidth[i];
        tmp1+=absNdot_nu[i]*sigma_HeI_nu(eVarrayCen[i])*eVarrayWidth[i];
        tmp2+=absNdot_nu[i]*sigma_HeII_nu(eVarrayCen[i])*eVarrayWidth[i];
    }
    time[0]=4*M_PI*pow(dist_cm,2)/tmp0;
    time[1]=4*M_PI*pow(dist_cm,2)/tmp1;
    time[2]=4*M_PI*pow(dist_cm,2)/tmp2;
    return 0;
}


int init_photo_bkg(cellProp *cell, double bkg[6]){

    if (cell->xHI<1e-15){cell->xHI=1e-15;}
    if (cell->xHeI<1e-15){cell->xHeI=1e-15;}
    if (cell->xHeII<1e-15){cell->xHeII=1e-15;}

    bkg[0] = recomb_rate(cell->T, HI)*cell->xHII*cell->ne/cell->xHI-cell->ne*eGamma_HI(cell->T);
    bkg[0] = bkg[0]>0 ? bkg[0] : 0;

    bkg[1] = recomb_rate(cell->T, HeI)*cell->xHeII*cell->ne/cell->xHeI-cell->ne*eGamma_HeI(cell->T);
    bkg[1] = bkg[1]>0 ? bkg[1] : 0;

    bkg[2] = recomb_rate(cell->T, HeII)*cell->xHeIII*cell->ne/cell->xHeII-cell->ne*eGamma_HeII(cell->T);
    bkg[2] = bkg[2]>0 ? bkg[2] : 0;
    
    bkg[3] = 0;
    bkg[4] = 0;
    bkg[5] = 0;
    return 0;
}
