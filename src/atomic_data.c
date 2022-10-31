#include<math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include"setup_path.h"
#include"atomic_data.h"
#include"base.h"

#define NUM_T_SAMPLE 1000
#define NUM_Z_SAMPLE 10001

//char dataDir[]="/Users/hqchen/reionization/pyfrit/pyfrit/RT1D/data/";
//extern char dataDir[];

double sigma_HI_nu(double nu){
    double nuHI=13.6;
    if (nu<nuHI)
        return 0;
    else
        return 6.30e-18*(1.34*pow((nu/nuHI),-2.99)-0.34*pow((nu/nuHI),-3.99));
}

double sigma_HeI_nu(double nu){
    double nuHeI=24.6;
    if (nu<nuHeI)
        return 0;
    else
        return 7.03e-18*(1.66*pow((nu/nuHeI),-2.05)-0.66*pow((nu/nuHeI),-3.05));
}

double sigma_HeII_nu(double nu){
    double nuHeII=54.4;
    if (nu<nuHeII)
        return 0;
    else
        return 1.50e-18*(1.34*pow((nu/nuHeII),-2.99)-0.34*pow((nu/nuHeII),-3.99));
}

double ion_potential(ionName species){
    switch(species){
        case HI:
            return 13.6;
        case HeI:
            return 24.6;
        case HeII:
            return 54.4;
        default:
            printf("ERROR: ion potential only for HI, HeI, HeII (not for HeI_r etc.)!!!\n");
            exit(VALUE_ERROR);

            
    }
} 

// ####  collisional-ionization rates [cm**3 s**-1] (names begin with eGamma)
double eGamma_HI(double T){
    if (T<10) {return 0;}
    double T5, eGamma;
    T5=T/1e5;
    eGamma=1.17e-10*pow(T,0.5)*exp(-157809.1/T)/(1+pow(T5,0.5));
    return eGamma;
}

double eGamma_HeI(double T){
    if (T<10) {return 0;}
    double T5, eGamma;
    T5=T/1e5;
    eGamma=4.76e-11*pow(T,0.5)*exp(-285335.4/T)/(1+pow(T5,0.5));
    return eGamma;
}

double eGamma_HeII(double T){
    if (T<10) {return 0;}
    double T5, eGamma;
    T5=T/1e5;
    eGamma=1.14e-11*pow(T,0.5)*exp(-631515.0/T)/(1+pow(T5,0.5));
    return eGamma;
}

double colli_ioniz_rate(double T, ionName species){
    switch(species){
        case HI:
            return eGamma_HI(T);
        case HeI:
            return eGamma_HeI(T);
        case HeII:
            return eGamma_HeII(T);
        default:
            printf("ERROR: collitional ionization rate only for HI, HeI, HeII (not for HeI_r etc.)!!!\n");
            exit(VALUE_ERROR);
    }
}

// ----- END --- collisional-ionization rates [cm**3 s**-1] -----

// ######## recombination rates [cm**3 s**-1] (names begin with alpha) #########

double recomb_rate_alphaHI_Abel97(double T){
    if (T<10) {T=10;}
    double Tbar=k_Boltzmann*T/eV2erg;
    double alpha= (exp(-28.6130338-0.72411256*log(Tbar)
            -2.02604473e-2*pow(log(Tbar),2) - 2.38086188e-3*pow(log(Tbar),3)
            -3.21260521e-4*pow(log(Tbar),4) - 1.42150291e-5*pow(log(Tbar),5)
            +4.98910892e-6*pow(log(Tbar),6) + 5.75561414e-7*pow(log(Tbar),7)
            -1.85676704e-8*pow(log(Tbar),8) - 3.07113524e-9*pow(log(Tbar),9)) );
    //printf("alpha HI=%e\n", alpha);
    return alpha;
}

double recomb_rate_alphaHeII_Abel97(double T){
    return (2* recomb_rate_alphaHI_Abel97(T/4));
}

double recomb_rate_alphaHeI_r_Abel97(double T){
    if (T<10) {T=10;}
    double Tbar=k_Boltzmann*T/eV2erg;
    double alpha= (3.925e-13*pow(Tbar,-0.6353));
    //printf("alphaHeI_r=%e\n", alpha);
    return alpha;
}

double recomb_rate_alphaHeI_d_Abel97(double T){
    if (T<1e3) {return 0;}
    double Tbar=k_Boltzmann*T/eV2erg;
    //printf("T=%e, Tbar=%e\n", T, Tbar);
    double alpha= ( 1.544e-9*pow(Tbar,-1.5)*exp(-48.596/Tbar)*(0.3+exp(8.1/Tbar)) );
    //printf("alphaHeI_d=%e\n", alpha);
    return alpha;
}
double recomb_rate_alphaHeI_Abel97(double T){
    return (recomb_rate_alphaHeI_r_Abel97(T)+recomb_rate_alphaHeI_d_Abel97(T));
}


gsl_spline *alphaA_HI_spline;
gsl_spline *alpha_r_HeI_spline;
gsl_spline *alpha_d_HeI_spline;
gsl_spline *alpha_HeI_spline;
gsl_spline *alpha_HeII_spline;
gsl_interp_accel *acc;

int init_load_recomb_data(){
    extern gsl_spline *alphaA_HI_spline;
    extern gsl_spline *alpha_r_HeI_spline;
    extern gsl_spline *alpha_d_HeI_spline;
    extern gsl_spline *alpha_HeI_spline;
    extern gsl_spline *alpha_HeII_spline;

    extern gsl_interp_accel *acc; 
    acc = gsl_interp_accel_alloc ();

    int nh;
    char fname[200];
    strcpy(fname,dataDir);
    strcat (fname,"recomb_rate_A_HI.dat");
    //printf("file name=%s",fname);
    int cnt=count_nrow_data_file(fname,'#',&nh);
    //printf("nrow=%d\n",cnt);
    
    double x[cnt], y[cnt];
    strcpy(fname,dataDir);
    strcat (fname,"recomb_rate_A_HI.dat");
    read_2D_data_file (fname, cnt, x, y, '#',nh);
    alphaA_HI_spline = gsl_spline_alloc (gsl_interp_cspline, NUM_T_SAMPLE);
    gsl_spline_init (alphaA_HI_spline, x, y, NUM_T_SAMPLE);

    strcpy(fname,dataDir);
    strcat (fname,"recomb_rate_r_HeI.dat");
    read_2D_data_file (fname, cnt, x, y, '#',nh);
    alpha_r_HeI_spline = gsl_spline_alloc (gsl_interp_cspline, NUM_T_SAMPLE);
    gsl_spline_init (alpha_r_HeI_spline, x, y, NUM_T_SAMPLE);

    strcpy(fname,dataDir);
    strcat (fname,"recomb_rate_d_HeI.dat");
    read_2D_data_file (fname, cnt, x, y, '#',nh);
    alpha_d_HeI_spline = gsl_spline_alloc (gsl_interp_cspline, NUM_T_SAMPLE);
    gsl_spline_init (alpha_d_HeI_spline, x, y, NUM_T_SAMPLE);

    strcpy(fname,dataDir);
    strcat (fname,"recomb_rate_HeI.dat");
    read_2D_data_file (fname, cnt, x, y, '#',nh);
    alpha_HeI_spline = gsl_spline_alloc (gsl_interp_cspline, NUM_T_SAMPLE);
    gsl_spline_init (alpha_HeI_spline, x, y, NUM_T_SAMPLE);

    strcpy(fname,dataDir);
    strcat (fname,"recomb_rate_HeII.dat");
    read_2D_data_file (fname, cnt, x, y, '#',nh);
    alpha_HeII_spline = gsl_spline_alloc (gsl_interp_cspline, NUM_T_SAMPLE);
    gsl_spline_init (alpha_HeII_spline, x, y, NUM_T_SAMPLE);

    return 0;
}

int free_loaded_recomb_data(){
    gsl_spline_free (alphaA_HI_spline);
    gsl_spline_free (alpha_r_HeI_spline);
    gsl_spline_free (alpha_d_HeI_spline);
    gsl_spline_free (alpha_HeI_spline);
    gsl_spline_free (alpha_HeII_spline);
    gsl_interp_accel_free (acc);
    return 0;
}

double recomb_rate(double T, ionName species){
    if (T<10) {
        T=10;
    } else if (T>1e8) {
        T=1e8;
    } else {
    }
    switch(species){
        case HI:
            //return gsl_spline_eval (alphaA_HI_spline, T, acc);
            return  recomb_rate_alphaHI_Abel97(T);
        case HeI:
            //return gsl_spline_eval (alpha_HeI_spline, T, acc);
            return recomb_rate_alphaHeI_Abel97(T);
        case HeII:
            //return gsl_spline_eval (alpha_HeII_spline, T, acc);
            return recomb_rate_alphaHeII_Abel97(T);
        case HeI_r:
            //return gsl_spline_eval (alpha_r_HeI_spline, T, acc);
            return recomb_rate_alphaHeI_r_Abel97(T);
        case HeI_d:
            //return gsl_spline_eval (alpha_d_HeI_spline, T, acc);
            return  recomb_rate_alphaHeI_d_Abel97(T);
    }
}



// ----- END ---  recombination rates [cm^3 s^-1]

// ##### Cooling Rates ###
double recomb_cooling(double T, double ne, double nHII, double nHeII, double nHeIII){
    if (T<10) {T=10;}
    double Lrec_HI, Lrec_HeI, Lrec_HeII;
    //printf("ne=%e, nHII= %e, T=%e, recrate=%e\n", ne, nHII, T, recomb_rate(T,HI));
    Lrec_HI=ne*nHII*1.036e-16*T*recomb_rate(T,HI);
    //printf("Lrec_HI=%e\n",Lrec_HI); 
    Lrec_HeI=ne*nHeII*(1.036e-16*T*recomb_rate(T,HeI_r)+6.526e-11*recomb_rate(T,HeI_d));
    //printf("Lrec_HeI=%e\n",Lrec_HeI); 
    Lrec_HeII=ne*nHeIII*1.036e-16*T*recomb_rate(T,HeII);
    //printf("Lrec_HeII=%e\n",Lrec_HeII); 
    //printf("finish?\n"); 
    //printf("Lrec_HI=%e, Lrec_HeI=%e, Lrec_HeII=%e\n", Lrec_HI, Lrec_HeI, Lrec_HeII);
    return (Lrec_HI+Lrec_HeI+Lrec_HeII);
}

double colli_ioniz_cooling(double T, double ne, double nHI, double nHeI, double nHeII){
    if (T<10) {T=10;}
    double Le_HI, Le_HeI, Le_HeII;
    Le_HI=ne*nHI*2.18e-11*colli_ioniz_rate(T,HI);
    Le_HeI=ne*nHeI*3.94e-11*colli_ioniz_rate(T,HeI);
    Le_HeII=ne*nHeII*8.72e-11*colli_ioniz_rate(T,HeII);
    //printf("Le_HI=%e, Le_HeI=%e, Le_HeII=%e\n", Le_HI, Le_HeI, Le_HeII);
    return (Le_HI+Le_HeI+Le_HeII);
}

double colli_excit_cooling(double T, double ne, double nHI, double nHeI, double nHeII){
    if (T<10) {return 0;}
    double T5, Lex_HI, Lex_HeI, Lex_HeII;
    T5=T/1e5;
    Lex_HI= ne*nHI*7.5e-19*exp(-118348/T)/(1+pow(T5,0.5));
    Lex_HeI=ne*ne*nHeI*9.10e-27*pow(T,-0.1687)*exp(-13179./T)/(1+pow(T5,0.5));
    Lex_HeII=ne*nHeII*5.54e-17*pow(T,-0.397)*exp(-473638)/(1+pow(T5,0.5));
    //printf("Lex_HI=%e, Lex_HeI=%e, Lex_HeII=%e\n", Lex_HI, Lex_HeI, Lex_HeII);
    return (Lex_HI+Lex_HeI+Lex_HeII);
}

double brem_cooling(double T, double ne, double nHII, double nHeII, double nHeIII){
    if (T<10) {T=10;}
    //double gff=1.5;
    double gff=1.1+0.34*pow(exp(-(5.5-log10(T))),2);
    double Lff;
    Lff=ne*(nHII+nHeII+4*nHeIII)*1.43e-27*pow(T,0.5)*gff;
    //printf("Lff=%e\n", Lff);
    return Lff;
}

double inv_comp_cooling(double T, double ne, double zuni){
    if (T<10) {T=10;}
    double Lc;
    Lc=5.65e-36*(T-2.73*(1+zuni))*pow((1+zuni),4)*ne;
    //printf("Lc=%e\n",Lc);
    return Lc;
}


double calc_Cooling(cellProp *cell, double zuni){
    double cooling_rate;
    double rec=recomb_cooling(cell->T, cell->ne, cell->nHII, cell->nHeII, cell->nHeIII);
    //printf("recomb_cooling=%e\n", rec);
    double colion=colli_ioniz_cooling(cell->T, cell->ne, cell->nHI, cell->nHeI, cell->nHeII);
    //printf("colli_ion_cooling=%e\n", colion);
    double colex=colli_excit_cooling(cell->T, cell->ne, cell->nHI, cell->nHeI, cell->nHeII);
    //printf("colli_ex_cooling=%e\n", colex);
    double brem=brem_cooling(cell->T, cell->ne, cell->nHII, cell->nHeII, cell->nHeIII);
    //printf("brem_cooling=%e\n", brem);
    double inv=inv_comp_cooling(cell->T, cell->ne, zuni);
    //printf("invC_cooling=%e\n", inv);
    //printf("rec=%e, colion=%e, colex=%e, brem=%e, inv=%e\n", rec, colion, colex, brem, inv);
    cooling_rate=rec+colion+colex+brem+inv;

    return cooling_rate;
}
